# detopt.py
# given a tree topology over CN-neutral mutations, assign CN-aberrant mutations 


### IMPORTS
import gurobipy as gp 
import os, re, math, copy, random, sys, json, time, argparse
import pandas as pd
from ast import literal_eval as lit
from util import *
from statistics import mean
from tree_scoring import anc_desc, diff_lin


### PARSER
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='DETOPT!', add_help=True)
	parser.add_argument('-d', '--dir',
						required=True,
						type=str,
						help='Data directory containing input files.')
	parser.add_argument('-i', '--input',
						required=True,
						type=str,
						help='Input .tree file with tumor subpopulation frequencies.')
	parser.add_argument('-r', '--reads',
						required=True,
						type=str,
						help='Input read data with subclonal copy number calls.')
	parser.add_argument('-g', '--ground_truth',
						required=False,
						type=str,
						help='Input ground truth tree file. To be used if running on simulations.')
	parser.add_argument('-o', '--out',
						required=False,
						type=str,
						help='Output file prefix. Defaults to prefix of the input file.')
	args = parser.parse_args()
	if not args.out:
		name = os.path.splitext(os.path.basename(args.input))[0]
		args.out = f'{name}.results'


### PARSE THE INPUTS
DATA_DIR = args.dir

print(f"INPUT FILE: {args.input}")

tree_file = f'{DATA_DIR}/{args.input}' #f'{TREE_DIR}/{SIM_NAME}.tree.tsv'.replace(f'h_{args.num_samples}', 'h_40')
true_tree_file = f'{DATA_DIR}/{args.ground_truth}'#f'{SIM_DIR}/{SIM_NAME}.tree.tsv'.replace('.neutral','').replace(f'h_{args.num_samples}', 'h_40')

tree = pd.read_csv(tree_file, header=0, index_col=False, sep='\t')
true_tree = pd.read_csv(true_tree_file, header=0, index_col=False, sep='\t')

reads_input = f'{DATA_DIR}/{args.reads}'
reads_df = pd.read_csv(reads_input, header=0, index_col=False, sep='\t')

C_MAX = 100 # large maximum copy number

# if on Biowulf or using Slurm job scheduler, try to grab number of CPUs
num_threads = os.environ.get('SLURM_CPUS_PER_TASK')
if not num_threads: 
	num_threads = 1
else: num_threads = int(num_threads)
print(f"Using {num_threads} threads.")


def freqs_from_true_tree(node, tree):
	res = {}
	samples = true_tree.iloc[0]['SAMPLE_IDS'].split(',')
	freqs = true_tree[true_tree['NODE_ID'] == node].squeeze()['SAMPLE_NODE_FREQUENCIES'].split(',')
	for i in range(len(samples)):
		res[samples[i]] = float(freqs[i])
	return res


def parse_tree(df, vafs, is_true_tree=False):
	samples = df.iloc[0]['SAMPLE_IDS'].split(',')
	vec, freqs = {}, {}
	for (i, row) in df.iterrows():
		node = row['NODE_ID']
		freqs[node] = {}
		vec[node] = row['PARENT_ID']
		nodefreqs = row['SAMPLE_NODE_FREQUENCIES'].split(',')
		samples = row['SAMPLE_IDS'].split(',')
		for i in range(len(samples)):
			sample = samples[i]
			freqs[node][sample] = float(nodefreqs[i])
	return vec, freqs, samples


def parse_data(df, mut):
	states, vafs = [], {}
	df = df[df['mut_index'] == mut]
	state_cols = [col for col in df.columns if 'state' in col and 'normal' not in col]
	for col in state_cols:
		state = tuple(df.iloc[0][col].split('|'))
		states.append(state)
	for (i, row) in df.iterrows():
		vaf = row['var_reads'] / (row['ref_reads'] + row['var_reads'])
		sample = row['sample']
		vafs[sample] = vaf
	return states, vafs


def fill(assignment, states, anc, nodes):
	'''
	given a single node assigned an aberrant CN state,
	fill in the copy numbers of all descendants of this node as the same state
	fill in all other nodes as the diploid state
	'''
	nodes_left = copy.deepcopy(nodes)
	nodes_left.remove(list(assignment.keys())[0])
	assert(len(assignment.keys()) == 1) # only root & another node assigned at the start
	assigned = list(assignment.keys())[0]
	assert(states[1] == assignment[assigned]) # hardcoded
	for node in descendants(assigned, anc):
		assignment[node] = assignment[assigned] # propagate the same state to descendants
		nodes_left.remove(node)
	for node in nodes_left: # remaining nodes get diploid state
		assignment[node] = states[0]
	assignment['ROOT'] = (1,1)
	return assignment


def solve(mutation, asg, parent, node_freqs, vafs):
	'''
	Input: mutation : mutation ID to attach to tree
		   asg : mapping of node -> copy number state tuple (A, B)
		   parent : dict mapping node -> parent
		   node_freqs : dict mapping node -> dict mapping sample -> node frequency
		   vafs : dict mapping sample id to observed vaf of mutation mutation
	Output:
		   score of a given assignment & gurobi model associated with it
	'''
	model = gp.Model("LP")
	model.Params.OutputFlag = 0
	model.Params.LogFile = ""
	model.Params.Threads = num_threads
	objective = 0.0

	delta = {}
	A = {} # number of copies of segment containing mut i in each node for allele A
	B = {} # number of copies of segment containing mut i in each node for allele A
	phi = {} # tells whether mutation is at A or B
	presence = {} # indicates whether mutation is present in a node
	U_A = {} # U_A[node] indicates if A[node] > 0 & all copies of A are mutated
	U_B = {} # U_B[node] indicates if B[node] > 0 & all copies of B are mutated

	nodes = parent.keys()
	for node in nodes: # delta_iv
		delta[node] = model.addVar(vtype=gp.GRB.BINARY, name=f'delta_{node}')
		A[node] = model.addVar(vtype=gp.GRB.INTEGER, name=f'A_{node}')
		B[node] = model.addVar(vtype=gp.GRB.INTEGER, name=f'B_{node}')
		presence[node] = model.addVar(vtype=gp.GRB.BINARY, name=f'pres_{node}')
		U_A[node] = model.addVar(vtype=gp.GRB.BINARY, name=f'U_A_{node}')
		U_B[node] = model.addVar(vtype=gp.GRB.BINARY, name=f'U_B_{node}')

	# there is never any mutational copy or presence at the ROOT node
	delta['ROOT'] = 0
	A['ROOT'] = 0
	B['ROOT'] = 0
	presence['ROOT'] = 0

	# define phi[allele] = 1 if mutation occurs at allele, 0 otherwise
	phi['A'] = model.addVar(vtype=gp.GRB.BINARY, name=f'phi_A')
	phi['B'] = model.addVar(vtype=gp.GRB.BINARY, name=f'phi_B')

	# CONSTRAINT 1 : mutation is assigned to a single node only
	model.addConstr(gp.quicksum([delta[node] for node in delta.keys()]) == 1)

	# CONSTRAINT 2 : mutation occurs either in allele A or B, but not both
	model.addConstr(phi['A'] + phi['B'] == 1)

	for node in nodes:
		if node == 'ROOT': continue

		# CONSTRAINTS 3+4 : mutated copies bounded by assignment & mutant allele
		model.addConstr(A[node] <= asg[node][0] * phi['A'])
		model.addConstr(B[node] <= asg[node][1] * phi['B'])

		# CONSTRAINTS 5+6 : mutational copies consistent with allelic & node assignment
		model.addConstr(A[node] >= delta[node] * phi['A'])
		model.addConstr(B[node] >= delta[node] * phi['B'])

		change_A = int(asg[node][0]) - int(asg[parent[node]][0])
		change_B = int(asg[node][1]) - int(asg[parent[node]][1])

		# CONSTRAINT 7+8: mutational copies are consistent if no CN event occurs
		if change_A == 0:
			model.addConstr(A[node] >= A[parent[node]])
			model.addConstr(A[node] <= A[parent[node]] + delta[node])		
		if change_B == 0:
			model.addConstr(B[node] >= B[parent[node]])
			model.addConstr(B[node] <= B[parent[node]] + delta[node])

		# CONSTRAINT 9 : mutational copies are consistent with loss events
		if change_A < 0:
			model.addConstr(A[node] <= A[parent[node]] + delta[node])
		if change_B < 0:
			model.addConstr(B[node] <= B[parent[node]] + delta[node])

		# CONSTRAINT 10 : mutational copies are consistent with gain events
		if change_A > 0:
			model.addConstr(A[node] <= A[parent[node]] + (change_A * presence[node]) + (change_A * delta[node]))
			model.addConstr(A[node] >= A[parent[node]])
			# CONSTRAINT 16 with U_A
			model.addConstr(A[node] >= A[parent[node]] + (change_A * U_A[node]))

		if change_B > 0:
			model.addConstr(B[node] <= B[parent[node]] + (change_B * presence[node]) + (change_B * delta[node]))
			model.addConstr(B[node] >= B[parent[node]])
			# CONSTRAINT 16 with U_B
			model.addConstr(B[node] >= B[parent[node]] + (change_B * U_B[node]))

		# CONSTRAINTS 11+12 : PRESENCE INDICATOR CONSTRAINTS
		model.addConstr(C_MAX * presence[node] >= A[node] + B[node])
		model.addConstr(presence[node] <= A[node] + B[node])

		# CONSTRAINTS 13+14+15 : if all allelic copies are mutated, any gain increases mutational copy number
		model.addConstr(U_A[node] >= A[node] - asg[node][0] + presence[node])
		model.addConstr(C_MAX * (1 - U_A[node]) >= asg[node][0] - A[node])
		model.addConstr(U_A[node] <= presence[node])

		model.addConstr(U_B[node] >= B[node] - asg[node][1] + presence[node])
		model.addConstr(C_MAX * (1 - U_B[node]) >= asg[node][1] - B[node])
		model.addConstr(U_B[node] <= presence[node])

	discrepancies = {}
	aux = {}
	for sample in vafs.keys():
		# define variant genomic content
		variant_genomic_content = gp.quicksum([presence[node]*node_freqs[node][sample]*(A[node]+B[node]) for node in nodes])
		# define total genomic content
		total_genomic_content = sum([(asg[node][0] + asg[node][1]) * node_freqs[node][sample] for node in nodes])

		# define auxiliary variables to allow absolute value in objective
		aux[sample] = model.addVar(vtype=gp.GRB.CONTINUOUS, name=f'aux_{sample}')

		discrepancy = (vafs[sample] - (variant_genomic_content / total_genomic_content))

		# enforce an absolute value objective
		model.addConstr(discrepancy <= aux[sample])
		model.addConstr(-discrepancy <= aux[sample])
		objective += aux[sample]

	model.setObjective(objective, gp.GRB.MINIMIZE)
	model.optimize()

	try:
		node_assigned = None
		for node in delta.keys():
			if node == 'ROOT': continue
			if delta[node].x == 1: node_assigned = node

		assert(node_assigned != None), f'Error: we assigned {mutation} to a node called None. Something is wrong with your .tree file.'
		return (node_assigned, model.objVal, model)
	except:
		# Model infeasible; failed to find any valid assignment - just return the root
		return ('ROOT', 1e6, model)



### Generate all possible CNA assignments, score each, and choose the best one
def assignCNAs_bruteforce(mutation, vec, node_freqs, states, vafs):
	'''
	generate all possible assignments of copy number states to tree nodes

	only the mutational (non-normal) states are given
	right now, we assume there are only two states

	score all possible assignments and choose the best one
	'''
	start = time.time()
	assignment_scores = {} # maps assigned node to its score from our qilp	
	anc = vec_to_anc_matrix(vec)
	nodes = list(vec.keys())
	states = states_to_int(states)
	for node in nodes:
		if node == 'ROOT': continue # don't assign mut to root node
		state_assignment = {node : states[1]} # pick the aberrant tumor state
		state_assignment = fill(state_assignment, states, anc, nodes) # fill in remaining states

		(node_assigned, objective, model) = solve(mutation, state_assignment, vec, node_freqs, vafs)

		hashed_assignment = json.dumps(state_assignment, sort_keys=True) # https://stackoverflow.com/questions/5884066
		assignment_scores[(node, node_assigned, hashed_assignment, model)] = objective

	'''
	for k in assignment_scores.keys():
		node_assigned = k[1]
		print(k[0], k[1], assignment_scores[k])
	'''
	(node_of_event, node_of_mutation, bestasg, bestmodel) = min(assignment_scores, key=assignment_scores.get)

	print(f"Optimal assignment: Mutation {mutation} attaches to Node {node_of_mutation}")
	end = time.time()
	bestasg = values_to_tuple(lit(bestasg)) # revert hashed assignment
	return node_of_event, node_of_mutation, bestasg, bestmodel



########
# MAIN #
########

muts_tested = sorted(sample_aberrant_muts(reads_df))
all_mutations = list(set(reads_df['mut_index'].values))

print("MUTS TESTED:", sorted(muts_tested))

if len(muts_tested) > 0: 
	correct_mut_assignments = 0
	avg_anc_desc = 0
	avg_diff_lin = 0
	ancdescs = []

	total_runtime = 0
	samples = tree.iloc[0]['SAMPLE_IDS'].split(',')
	vafs = parse_vafs(reads_df, samples)
	vec, freqs, samples_chosen = parse_tree(tree, vafs)
	true_vec, _, _ = parse_tree(true_tree, vafs, True)

	final_assignments = {}
	for mut in muts_tested:
		# parse the data and solve the model for each mutation individually
		states, vafs_of_mut = parse_data(reads_df, mut)
		node_of_event, node_of_mutation, bestasg, bestmodel = assignCNAs_bruteforce(mut, vec, freqs, states, vafs_of_mut)

		final_assignments[mut] = (node_of_mutation, node_of_event)

		if args.ground_truth is not None: # score simulations
			actual_node = get_actual_node(mut, tree, true_tree)
			lab = 'SNV_IDS'

			ancDesc = anc_desc(mut, node_of_mutation, all_mutations, vec_to_anc_matrix(true_vec),
							   vec_to_anc_matrix(vec), true_tree, tree, 'SNV_IDS')
			diffLin = diff_lin(mut, node_of_mutation, all_mutations, vec_to_anc_matrix(true_vec),
							   vec_to_anc_matrix(vec), true_tree, tree, 'SNV_IDS')
			print(f"ancDesc: {ancDesc}, diffLin: {diffLin}")
			avg_anc_desc += ancDesc
			ancdescs.append(ancDesc)
			avg_diff_lin += diffLin

	if args.ground_truth is not None:
		print(f"Ancestor-Descendant Accuracy: {avg_anc_desc / len(muts_tested)}")
		print(f"Different-Lineage Accuracy: {avg_diff_lin / len(muts_tested)}")

	with open(f'{args.out}.tsv', 'w') as o:
		o.write('MUT_ID\tNODE_ASSIGNED\tNODE_OF_CNA\n')
		for k in final_assignments.keys():
			node_of_mut, node_of_cna = final_assignments[k]
			o.write(f'{k}\t{node_of_mut}\t{node_of_cna}\n')

else:
	print("There were no aberrant mutations to place. ")