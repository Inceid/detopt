# detopt
`DETOPT` is a method to efficiently assign SNVs affected by copy number aberrations trees of tumor evolution.


## Dependencies
`DETOPT` Requires the following dependencies:
- Gurobi version 9.0.1 or later
- Python version 3.6 or later


## Setup
For this version, we assume the user is using Linux. Setting up `DETOPT` on your local machine can be accomplished with the following stepsL

### 1. Resolving Dependencies

- Gurobi offers free software licenses for academic use. Download Gurobi [here](https://www.gurobi.com/downloads/gurobi-optimizer-eula/) and follow the instructions [here](https://www.gurobi.com/academia/academic-program-and-licenses/) for obtaining a license.
- You can download and install Python [here](https://www.python.org/downloads/).


### 2. Download

Use this command to download `DETOPT`'s source code to your machine.

```
git clone https://github.com/Inceid/detopt
```

### 3. Testing

You can run `DETOPT` on the provided simulated dataset with the following command. Please ensure that any simulated `reads.input` files are contained in the `sims/` directory in order for this command to work.

```
python detopt.py -d sims -i simNo_1-s_21-m_100-h_10-cna_0.1.tree -r simNo_1-s_21-m_100-h_10-cna_0.1-cov_100.reads.input -g simNo_1-s_21-m_100-h_10-cna_0.1.groundtruth.tree
```

### 4. Input

`DETOPT` requires two inputs to run on a simulated patient file: a `.tree` file containing information from the dataset's phylogeny, and a `.reads.input` file containing observed read counts and copy number calls.

The `.tree` file describes the essential components of an inferred phylogeny for `DETOPT`'s functionality. It should contain columns specifying each tree node's `NODE_ID`, `PARENT_ID`, `SNV_IDS`, `SAMPLE_IDS`, and `SAMPLE_NODE_FREQUENCIES`. These are defined as follows.

- `NODE_ID`: a name for a mutational subclone. We require that the ID of the root node is 'ROOT'. 
- `PARENT_ID`: a name for the parent node of a given node. We require that the `ROOT` have no parents, or parent `NONE`.
- `SNV_IDS`: a name for the SNVs (mutations) assigned to a given node.
- `SAMPLE_IDS`: a comma-separated string of all sample IDs.
- `SAMPLE_NODE_FREQUENCIES`: a comma-separated string of the subpopulation frequency of a node in each sample, listed in the same order as the sample IDs in `SAMPLE_IDS`. For example, expect the 1st value of `SAMPLE_NODE_FREQUENCIES` to correspond to the node's frequency in the 1st sample in `SAMPLE_IDS`.


The `.reads.input` file describes the essential components of mutational and copy number data obtained from bulk whole-exome/genome sequencing following by mutation and copy number calling. We use an in-house pipeline for joint multi-sample variant calling and `HATCHet` for joint multi-sample copy number calling.

Upon obtaining VCFs and HATCHet's `seg.ucn` output files, the user should convert the data to the following format (we only show one SNV, `mut0`).

	mut_index	sample	var_reads	ref_reads	normal_state	normal_prop	tumor1_state	tumor1_prop	tumor2_state	tumor2_prop
	mut0	S0	20	30	1|1	0.1	1|1	0.2	2|1	0.7
	mut0	S1	10	40	1|1	0.1	1|1	0.4	2|1	0.5
	mut0	S2	0	50	1|1	1.0	1|1	0.0	2|1	0.0


Please see the `sims/` directory for an example of both inputs. The user should provide read counts for each mutation for each sample, even if the variant read counts are 0 in some samples.

The user should be able to run `DETOPT` on a real dataset with the command line options `-d (--dir)`, `-i (--input)`, and `-r (--reads)` specifying their corresponding inputs properly. Option `-g (--ground_truth)` should only be used when testing simulations.


### 5. Output

By default, `DETOPT` outputs a single file with suffix `results.tsv` which, for each mutation afflicted by copy number aberrations, reports the mutation ID, the phylogenetic node to which it was assigned, and the phylogenetic node to which its associated aberrant state(s) were assigned. An example of this format is given below.

	MUT_ID	NODE_ASSIGNED	NODE_OF_CNA
	mut11	2	7
	mut2	15	0
	mut33	17	13
	mut38	6	10
	mut50	17	15
	mut61	11	6
	mut69	10	15
	mut78	3	1
	mut85	16	19
	mut86	16	0


### 6. Correspondence

Feel free to open an issue if you'd like to discuss usage of `DETOPT`.
