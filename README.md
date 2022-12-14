# detopt
`DETOPT` is a method to efficiently assign SNVs affected by copy number aberrations to trees of tumor progression (also known as trees of tumor evolution).


## Dependencies
`DETOPT` Requires the following dependencies:
- Gurobi version 9.0.1 or later
- Python version 3.6 or later


## Setup
`DETOPT` is currently supported on LINUX OS. Setting up `DETOPT` on your local machine can be accomplished with the following steps:

### 1. Resolving Dependencies

- Gurobi offers free software licenses for academic use. Download Gurobi here: www.gurobi.com/downloads/gurobi-optimizer-eula/. 
- Follow the instructions at www.gurobi.com/academia/academic-program-and-licenses/ for obtaining a license.
- Note that successful use of Gurobi requires modification of multiple of environment variables as described in the Software Installation Guide here: https://www.gurobi.com/documentation/9.5/quickstart_linux/software_installation_guid.html.
- For downloading and installing Python please visit: https://www.python.org/downloads/.
- Finally, download and install the Gurobi Python interface `gurobipy` with:

```
python -m pip install gurobipy
```

### 2. Download

Use this command to download `DETOPT`'s source code to your machine.

```
git clone https://github.com/Inceid/detopt
```

### 3. Testing

You can run `DETOPT` on the provided simulated dataset with the following command. Please ensure that any simulated `reads.input` files are contained in the `input/` directory in order for this command to work.

```
python detopt.py -d input -i example.tree -r example.reads.input
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


Please see the `input/` directory for an example of both inputs. The user should provide read counts for each mutation for each sample, even if the variant read counts are 0 in some samples.

The user should be able to run `DETOPT` on a real dataset with the command line options `-d (--dir)`, `-i (--input)`, and `-r (--reads)` specifying their corresponding inputs properly. Option `-g (--ground_truth)` should only be used when testing simulations.


### 5. Output

By default, `DETOPT` outputs a single file with suffix `results.tsv` which, for each mutation afflicted by copy number aberrations, reports the mutation ID, the phylogenetic node to which it was assigned, and the phylogenetic nodes to which each of its associated aberrant state(s) were assigned. An example of this format is given below. The `CNA_ASSIGNMENTS` column can be interpreted as follows (for example): `0:(1,1);1:(2,1)` indicates Node 0 maps to copy number state (1,1), whereas Node 1 maps to copy number state (2,1).

	MUT_ID	NODE_ASSIGNED	CNA_ASSIGNMENTS
	mut11	2	0:(1,1);1:(1,1);10:(1,0);11:(1,0);12:(1,0);13:(1,1);14:(1,1);15:(1,1);16:(1,0);17:(1,1);18:(1,1);19:(1,1);2:(1,1);3:(1,1);4:(1,1);5:(1,1);6:(1,1);7:(1,0);8:(1,0);9:(1,0);ROOT:(1,1)
	mut2	15	0:(2,1);1:(2,1);10:(2,1);11:(2,1);12:(2,1);13:(2,1);14:(2,1);15:(2,1);16:(2,1);17:(2,1);18:(2,1);19:(2,1);2:(2,1);3:(2,1);4:(2,1);5:(2,1);6:(2,1);7:(2,1);8:(2,1);9:(2,1);ROOT:(1,1)
	mut33	17	0:(1,1);1:(1,1);10:(1,1);11:(1,1);12:(1,1);13:(2,1);14:(1,1);15:(1,1);16:(1,1);17:(1,1);18:(1,1);19:(1,1);2:(1,1);3:(1,1);4:(1,1);5:(1,1);6:(1,1);7:(1,1);8:(1,1);9:(1,1);ROOT:(1,1)
	mut38	6	0:(1,1);1:(1,1);10:(1,0);11:(1,1);12:(1,1);13:(1,1);14:(1,1);15:(1,1);16:(1,1);17:(1,1);18:(1,1);19:(1,1);2:(1,1);3:(1,1);4:(1,1);5:(1,1);6:(1,1);7:(1,1);8:(1,1);9:(1,1);ROOT:(1,1)
	mut50	17	0:(1,1);1:(1,1);10:(1,1);11:(1,1);12:(1,1);13:(1,1);14:(1,1);15:(2,1);16:(1,1);17:(1,1);18:(1,1);19:(1,1);2:(1,1);3:(1,1);4:(1,1);5:(1,1);6:(1,1);7:(1,1);8:(1,1);9:(1,1);ROOT:(1,1)
	mut61	11	0:(1,1);1:(1,1);10:(1,0);11:(1,0);12:(1,0);13:(1,1);14:(1,1);15:(1,1);16:(1,0);17:(1,1);18:(1,1);19:(1,0);2:(1,1);3:(1,1);4:(1,1);5:(1,1);6:(1,0);7:(1,0);8:(1,0);9:(1,0);ROOT:(1,1)
	mut69	10	0:(1,1);1:(1,1);10:(1,1);11:(1,1);12:(1,1);13:(1,1);14:(1,1);15:(2,1);16:(1,1);17:(1,1);18:(1,1);19:(1,1);2:(1,1);3:(1,1);4:(1,1);5:(1,1);6:(1,1);7:(1,1);8:(1,1);9:(1,1);ROOT:(1,1)
	mut78	3	0:(1,1);1:(2,1);10:(2,1);11:(2,1);12:(2,1);13:(2,1);14:(1,1);15:(1,1);16:(2,1);17:(1,1);18:(1,1);19:(2,1);2:(1,1);3:(1,1);4:(1,1);5:(1,1);6:(2,1);7:(2,1);8:(2,1);9:(2,1);ROOT:(1,1)
	mut85	16	0:(1,1);1:(1,1);10:(1,1);11:(1,1);12:(1,1);13:(1,1);14:(1,1);15:(1,1);16:(1,1);17:(1,1);18:(1,1);19:(2,1);2:(1,1);3:(1,1);4:(1,1);5:(1,1);6:(1,1);7:(1,1);8:(1,1);9:(1,1);ROOT:(1,1)
	mut86	16	0:(2,1);1:(2,1);10:(2,1);11:(2,1);12:(2,1);13:(2,1);14:(2,1);15:(2,1);16:(2,1);17:(2,1);18:(2,1);19:(2,1);2:(2,1);3:(2,1);4:(2,1);5:(2,1);6:(2,1);7:(2,1);8:(2,1);9:(2,1);ROOT:(1,1)



### 6. Correspondence

If you encounter any issues during setting up or running `DETOPT`, please open an issue on the repository or contact Suraj Joshi (`suraj.joshi@nih.gov`) or Salem Malikic (`salem.malikic@nih.gov`).
