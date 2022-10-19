# ancestor-descendant scoring

def splitmuts(muts):
	if isinstance(muts, float): return [] # for NaN cases
	return muts.split(',')


def find_node(mut, tree, mut_col_name):
	tree2 = tree.copy(deep=True)
	tree2[mut_col_name] = tree2[mut_col_name].map(splitmuts)
	node = 'ROOT'
	for (i, row) in tree2.iterrows():
		if mut in row[mut_col_name]:
			node = row['NODE_ID']
	return node


def anc_desc(aberrant_mut, aberrant_mut_node,
			 muts, anc_ground_truth, anc_inferred, 
			 ground_truth, inferred, mut_col_name):
	'''
	input: aberrant mutation ID
		   all mutation IDs
		   ancestry matrix of ground truth tree
		   ancestry matrix of inferred tree
		   ground truth tree df
		   inferred tree df
	output: ancestor-descendant accuracy score.
	'''
	print(f"Aberrant mut: {aberrant_mut}, Aberrant node: {aberrant_mut_node}")
	total = 0
	success = 0
	for mut in muts: # over all mutations
		gt_node_i = find_node(aberrant_mut, ground_truth, 'SNV_IDS')
		gt_node_j = find_node(mut, ground_truth, 'SNV_IDS')

		inf_node_i = aberrant_mut_node
		inf_node_j = find_node(mut, inferred, mut_col_name)

		if gt_node_i == gt_node_j:
			continue

		elif (anc_ground_truth[gt_node_j][gt_node_i] == 0 and 
			  anc_ground_truth[gt_node_i][gt_node_j] == 0):
			continue

		elif anc_ground_truth[gt_node_j][gt_node_i]:
			total += 1
			if anc_inferred[inf_node_j][inf_node_i]:
				success += 1

		elif anc_ground_truth[gt_node_i][gt_node_j]:
			total += 1
			if anc_inferred[inf_node_i][inf_node_j]:
				success += 1

		else:
			assert False, "Error! There's some mysterious bug."
	return success / total


def diff_lin(aberrant_mut, aberrant_mut_node,
			 muts, anc_ground_truth, anc_inferred, 
			 ground_truth, inferred, mut_col_name):
	'''
	input: aberrant mutation ID
		   all mutation IDs
		   ancestry matrix of ground truth tree
		   ancestry matrix of inferred tree
		   ground truth tree df
		   inferred tree df
	output: different-lineage accuracy score.
	'''
	total = 0
	success = 0
	for mut in muts:
		gt_node_i = find_node(aberrant_mut, ground_truth, 'SNV_IDS')
		gt_node_j = find_node(mut, ground_truth, 'SNV_IDS')

		inf_node_i = aberrant_mut_node
		inf_node_j = find_node(mut, inferred, mut_col_name)
		
		if gt_node_i == gt_node_j:
			continue

		elif (anc_ground_truth[gt_node_j][gt_node_i] == 0 and 
			  anc_ground_truth[gt_node_i][gt_node_j] == 0):
			total += 1
			if (anc_inferred[inf_node_i][inf_node_j] == 0 and
				anc_inferred[inf_node_j][inf_node_i] == 0):
				success += 1

		elif anc_ground_truth[gt_node_j][gt_node_i]:
			continue

		elif anc_ground_truth[gt_node_i][gt_node_j]:
			continue

		else:
			assert False, "Error! There's some mysterious bug."

	if total == 0: return 1.0
	return success / total




'''
take copy number affected variant i.
total = 0
success = 0 
For each other variant j:
       if, in the ground truth, j is at same node or at the node that is on different branch compared to i:
             continue
       elif, in the ground truth, j is ancestor of i:
                    total = total + 1
                    if, in the inferred tree, j is ancestor of i: 
                             success = success + 1
       elif, in the ground truth tree, i is ancestor of j:
                    total = total + 1
                    if, in the inferred tree, i is ancestor of j:
                              success = success + 1
       else: 
                   assert False, "Oooops. There is some bug in Salem's pseudocode or in my code"


accuracy(i) = success/total
(total should always be non-zero but if you encounter division by zero error)

Do the same for different lineage, this time only considering variants that are on different lineages (different tree branches).

To get these fast, I recommend you to first caluclate anc_matrix[node_1][node_2]  for the ground truth tree and for the inferred tree. Then from there you can quickly compute these. Just be careful when treating variants from the same node to not introduce bug as that is edge case.
'''