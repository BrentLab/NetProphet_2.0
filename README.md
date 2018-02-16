# NETPROPHET 2.0

NetProphet 2.0 is a second-generation “data light” TF-network mapping algorithm. It requires only data that can be generated from low-cost, reliable, and easily scalable experimental methods. NetProphet 2.0 relies on three fundamental ideas. First, combining several expression-based network algorithms that use different types of models can yield better results than using either one alone. Second, TFs with similar DNA binding domains (in terms of amino acid sequence) tend to bind similar sets of target genes. Third, even an imperfect net-work map can be used to infer models of each TF’s DNA binding prefer-ences from the promoter sequences of its putative targets and these mod-els can be used to further refine the network.

![NetProphet2.0_overview](NP2_overview.png)

### SYSTEM REQUIREMENTS

* Slurm workload manager (tested on v15.08.7)
* Open MPI (tested on v1.8.8)
* R (>= v3.2, tested on v3.2.1)
* Python (>= v2.7, tested on v2.7.10)
* Python (>= v3.4, tested on 3.4.3+)

### INSTALLATION INSTRUCTIONS

1. Configure NetProphet 2.0 directory
	
	```
	export NETPROPHET2_DIR=/path/to/NetProphet_2.0
	```

2. Install Snakemake (workflow management system)

	```
	cd ${NETPROPHET2_DIR}/SRC/
	tar -zxvf snakemake-3.8.2.tar.gz
	cd snakemake-3.8.2/
	python3 setup.py build
	python3 setup.py install --user
	```

3. Install FIRE program

	```
	cd ${NETPROPHET2_DIR}/SRC/
	unzip -q FIRE-1.1a.zip
	cd FIRE-1.1a/
	chmod 775 configure
	./configure
	make
	```

4. Install MEME suite

	```	
	cd ${NETPROPHET2_DIR}/SRC/
    mkdir -p meme/
	tar -zxvf meme_4.9.1.tar.gz
	cd meme_4.9.1/
	./configure --prefix=${NETPROPHET2_DIR}/SRC/meme --with-url="http://meme.nbcr.net/meme"
	make
	make test
	make install
	```

5. Install R packages

	```
	cd ${NETPROPHET2_DIR}/SRC/R_pkgs/
	module load R/3.2.1  #if not loaded by default
	```
	
	1. Install the following package version: 
		* lars v0.9-8
		* BayesTree v0.3-1.3

	```R
	> install.packages("BayesTree_0.3-1.3.tar.gz", lib="<your_local_R_lib>") 
	> install.packages("lars_0.9-8.tar.gz", lib="<your_local_R_lib>")
	```
	
	2. **[Optional]** Install Rmpi, if not available in your system.
	```R
	> install.packages("Rmpi_0.5-9.tar.gz", lib="<your_local_R_lib>")
	```

### EXAMPLE USAGE

1. Set pipeline paths. Execute the following for each NetProphet 2.0 run, or add them to `$HOME/.profile` or `$HOME/.bash_profile` once for all.
	
	```".bash_profile"
	export R_LIBS=<your_local_R_lib>
	export PATH=$HOME/.local/bin:$PATH
	export NETPROPHET2_DIR=/path/to/NetProphet_2.0
	export PATH=${NETPROPHET2_DIR}:$PATH
	export FIREDIR=${NETPROPHET2_DIR}/SRC/FIRE-1.1a/
	export PATH=${FIREDIR}:$PATH
	export PATH=${NETPROPHET2_DIR}/SRC/meme/bin/:$PATH
	```

2. Configure input, parameter and output in `config.json`. All required resources files are in directory `RESOURCES/`, and the output network file is in directory `OUTPUT/`. See **DESCRIPTION OF RESOURCE FILES**, **DESCRIPTION OF INPUT PARAMETERS** and **DESCRIPTION OF OUTPUT FILE** for details.
	
	```json
	{
		"NETPROPHET2_DIR": "/path/to/NetProphet_2.0",
		"FILENAME_EXPRESSION_DATA": "data.expr",
		"FILENAME_DE_ADJMTR": "signed.de.adj",
		"FILENAME_GENES": "genes",
		"FILENAME_REGULATORS": "regulators",
		"FILENAME_SAMPLE_CONDITIONS": "conditions",
		"DIR_DBD_PID": "DBD_PIDS",
		"FILENAME_PROMOTERS": "promoter.fasta",
		"PROMOTER_LENGTH": 600,
		"MOTIF_THRESHOLD": 16,
		"FILENAME_NETPROPHET2_NETWORK": "netprophet2_network.adjmtr"
	}
	```

3. Run NetProphet 2.0.

	```
	cd ${NETPROPHET2_DIR}
	sbatch NetProphet2
	```

4. Progress monitoring and debugging. After the SLURM job starts, a log file records the progress as the following.
	
	```
	Unlocking working directory.
	Provided cores: 2
	Rules claiming more threads will be scaled down.
	Job counts:
		count	jobs
		1	all
		1	assemble_final_network
		1	build_motif_network
		1	combine_npwa_bnwa
		1	infer_motifs
		1	make_directories
		1	map_bart_network
		1	map_np_network
		1	prepare_resources
		1	score_motifs
		1	weighted_average_bart_network
		1	weighted_average_np_network
		12
	rule make_directories:
		output : ...
	1 of 12 steps (8%) done
	rule prepare_resources:
	...
	2 of 12 steps (17%) done
	...
	...
	...
	12 of 12 steps (100%) done
	```

> NOTE: The resource data provided is used for mapping a Yeast subnetwork. Visit http://mblab.wustl.edu/software.html for the resources for mapping whole TF network in yeast and fruit fly.

### DESCRIPTION OF RESOURCE FILES

FILE/DIRECTORY | DESCRITPION
--- | ---
FILENAME_EXPRESSION_DATA | A matrix of the log2 fold-change expression values in the samples with respect to those in wildtype. Rows represent genes, columns represent samples/conditions, i.e. the matrix dimension is number of genes x number of samples.
FILENAME_DE_ADJMTR | A adjacency matrix of the interactions between regualtors and target genes, which are calculated via differential expression analysis. The rows represent regulators/TFs and the columns represent genes, i.e. the matrix dimension is number of regulators x number of target genes. For each possible interaction between regulator i (Ri) and target gene j (Tj), set entry Mij to the signed logged differential expression significance of Tj when Ri is perturbed. If Ri has not been perturbed, then set Mij = 0 for all j. See **CALCULATING THE DIFFERENTIAL EXPRESSION COMPONENT** for details.
FILENAME_GENES | A list of gene names. Capitalized systematic names are recommended.
FILENAME_REGULATORS | A list of gene names that encode transcription factors (TFs). These regulators must be included in the list of gene names. The regulator names should have the same naming scheme as the gene names. 
FILENAME_SAMPLE_CONDITIONS | A list of samples/conditions. If a gene was perturbed in a condition, set the condition name as the gene name; otherwise, set as any identifier without space delimiter.
FILENAME_PROMOTERS | The promoter sequences of the target genes in Fasta format. The header of each promter is the gene name only.
DIR_DBD_PID | A directory of the percent identities (PIDs) between the DNA binding domains (DBDs). Each file is titled as the name of the regulator associated with a DBD. There are two columns in the file: each entry of the first column is the regulator name associated with other DBDs, and the entry of the second column is the corresponding PID calculated beforehand. See **CALCULATING THE PERCENT IDENTITIES BETWEEN THE DBDS** for details.

### DESCRIPTION OF INPUT PARAMETERS

PARAMETER | DESCRITPION
--- | ---
PROMOTER_LENGTH | The number of base paires of the promoter used.
MOTIF_THRESHOLD | The threshold on the robustness score calculated in FIRE's jack-knife valdiation. Choose a value between 16 and 20. NetProphet 2.0 paper used the threshold of 20.  

### DESCRIPTION OF OUTPUT FILE

FILE | DESCRITPION
--- | ---
netprophet2_network.adjmtr | A adjacency matrix of the final scores predicted by NetProphet 2.0. The rows represent regulators/TFs and the columns represent genes, i.e. the matrix dimension is # of regulators x # of target genes. Each entry Mij of matrix M is the score of the interaction between regulator Ri and target gene Tj. In this matrix, interactions with higher scores are more likely to be direct, functional interactions.

### CALCULATING THE DIFFERENTIAL EXPRESSION COMPONENT

##### Microarray expression profiling data

For each TF perturbation, for each gene in the perturbation condition, we recommend that you use LIMMA to calculate the log odds that the gene is differentially expressed in the perturbation condition compared to the wild type (WT) condition. The differential expression component is a signed confidence score Dij, which is calculated using the log odds score Li(j) and the log2-fold change Yi(j) of gene j and TF i as follows.

Dij =  Li(j)*sgn(Yi(j) when Li(j) > 0 and Dij =  0 when Li(j) <= 0

##### RNA-Seq expression profiling data

For each TF perturbation, we recommend that you use Cuffdiff to calculate the significance of differential expression (i.e. the uncorrected p-value and the FDR-adjusted p-value) of each gene in the perturbation condition compared to the WT condition. The differential expression component is a signed confidence score Dij, which is calculated using the uncorrected p-value Pi(j), the FDR-adjusted p-value Fi(j), and the log2-fold change Yi(j) of gene j and TF i as follows.
	
Dij =  -ln(Pi(j))*sgn(Yi(j) when Fi(j) <= 0.05 and Dij =  0 when Fi(j) >= 0.05

### CALCULATING THE PERCENT IDENTITIES BETWEEN THE DBDS

See supplemental package at https://github.com/yiming-kang/DBD_PercentIdentity_Calculation for details.

### REFERENCES

Kang, Y, et al. NetProphet 2.0: Mapping Transcription Factor Networks by Exploiting Scalable Data Resources. Bioinformatics 2018;34(2):249–257.

Haynes, B.C., et al. Mapping functional transcription factor networks from gene expression data. Genome research 2013;23(8):1319-1328.

Chipman, H.A., George, E.I. and McCulloch, R.E. BART: Bayesian additive regression trees. 2010:266-298.

Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.

Elemento, O., Slonim, N. and Tavazoie, S. A universal framework for regulatory element discovery across all genomes and data types. Mol Cell 2007;28(2):337-350.

Grant, C.E., Bailey, T.L. and Noble, W.S. FIMO: scanning for occurrences of a given motif. Bioinformatics (Oxford, England) 2011;27(7):1017-1018.

Smyth GK. 2005. "Limma: Linear models for microarray data". Bioinformatics and computational biology solutions using R and Bioconductor (ed. Gentleman R, et al.), pp. 397–420. Springer, New York.

Efron, Bradley; Hastie, Trevor; Johnstone, Iain; Tibshirani, Robert. Least angle regression. Ann. Statist. 32 (2004), no. 2, 407--499.

Yu, Hao. "Rmpi: parallel statistical computing in R." R News 2.2 (2002): 10-14.

