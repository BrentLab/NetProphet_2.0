workdir: config["NETPROPHET2_DIR"]

NETPROPHET2_NETWORK = "/".join([config["OUTPUT_DIR"],config["FILENAME_NETPROPHET2_NETWORK"]])

rule all:
	input:
		NETPROPHET2_NETWORK

rule make_directories:
	output:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],"tmp"]),
		n = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],"networks/"]),
		m = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],"motif_inference/"]),
		s = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"motif_inference/network_scores/"]),
		b = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"motif_inference/network_bins/"]),
		p = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"motif_inference/motifs_pfm/"]),
		q = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"motif_inference/motifs_score/"]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.make_dir"])
	shell:
		"""
		printf "\e[1mStep 0: Initializing\e[0m\n"; mkdir -p {output.r}; mkdir -p {output.n}; mkdir -p {output.m}; mkdir -p {output.s}; mkdir -p {output.b}; mkdir -p {output.p}; mkdir -p {output.q}; touch {output.flag}
		"""

rule prepare_resources:
	input:
		g = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_GENES"]]),
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		e = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_EXPRESSION_DATA"]]),
		c = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_SAMPLE_CONDITIONS"]]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.make_dir"])
	output:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/rdata.expr"]),
		f = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/data.fc.tsv"]),
		a = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/allowed.adj"]),
		p1 = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/data.pert.adj"]),
		p2 = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/data.pert.tsv"]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.prepare_resources"])
	shell:
		"""
		python CODE/prepare_resources.py -g {input.g} -r {input.r} -e {input.e} -c {input.c} -or {output.r} -of {output.f} -oa {output.a} -op1 {output.p1} -op2 {output.p2}; touch {output.flag}; printf "\e[1m[Step 0 completed]\e[0m\n\n";
		"""

rule map_np_network:
	input:
		g = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_GENES"]]),
		f = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		t = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_EXPRESSION_DATA"]]),
		d = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_DE_ADJMTR"]]),
		u = "/".join([config["NETPROPHET2_DIR"],"SRC/NetProphet1/"]),
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/rdata.expr"]),
		a = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/allowed.adj"]),
		p = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/data.pert.adj"]),
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/"]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.prepare_resources"])
	output:
		n = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/np.adjmtr"])
	shell:
		"""
		printf "\e[1mStep 1: Mapping NetProphet 1.0 network\e[0m\n"; ./SRC/NetProphet1/netprophet -m -u {input.u} -t {input.t} -r {input.r} -a {input.a} -p {input.p} -d {input.d} -g {input.g} -f {input.f} -o {input.o} -n {output.n}; printf "\e[1m[Step 1 completed]\e[0m\n\n";
		"""

rule map_bart_network:
	input:
		f = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		t = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/data.fc.tsv"]),
		p = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/data.pert.tsv"]),
		n = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/np.adjmtr"])
	output:
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/bn.adjmtr"])
	shell:
		"""
		printf "\e[1mStep 2: Mapping BART network\e[0m\n"; ./CODE/run_build_bart_network.sh {input.t} {input.p} {input.f} {output.o} true; printf "\e[1m[Step 2 completed]\e[0m\n\n";
		"""

rule weighted_average_np_network:
	input:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		a = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["DBD_PID_DIR"]]),
		n = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/np.adjmtr"]),
		b = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/bn.adjmtr"])
	output:
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/npwa.adjmtr"])
	shell:
		"""
		printf "\e[1mStep 3.1: Weighted averaging NetProphet 1.0 scores\e[0m\n"; python CODE/weighted_avg_similar_dbds.py -n {input.n} -r {input.r} -a {input.a} -d 50 -f single_dbds -o {output.o}; printf "\e[1m[Step 3.1 completed]\e[0m\n\n";
		"""

rule weighted_average_bart_network:
	input:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		a = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["DBD_PID_DIR"]]),
		n = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/bn.adjmtr"]),
		b = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/npwa.adjmtr"])
	output:
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/bnwa.adjmtr"])
	shell:
		"""
		printf "\e[1mStep 3.2: Weighted averaging BART scores\e[0m\n"; python CODE/weighted_avg_similar_dbds.py -n {input.n} -r {input.r} -a {input.a} -d 50 -f single_dbds -o {output.o}; printf "\e[1m[Step 3.2 completed]\e[0m\n\n";
		"""

rule combine_npwa_bnwa:
	input:
		n = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/npwa.adjmtr"]),
		b = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/bnwa.adjmtr"])
	output:
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/npwa_bnwa.adjmtr"])
	shell:
		"""
		printf "\e[1mStep 4: Combining networks\e[0m\n"; Rscript CODE/quantile_combine_networks.r {input.n} {input.b} {output.o}; printf "\e[1m[Step 4 completed]\e[0m\n\n";
		"""

rule infer_motifs:
	input:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		t = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_GENES"]]),
		p = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_PROMOTERS"]]),
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],""]),
		a = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/npwa_bnwa.adjmtr"])
	output:
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.infer_motifs"])
	shell:
		"""
		printf "\e[1mStep 5.1: Inferring TF binding motifs\e[0m\n"; ./CODE/run_infer_motifs.sh {input.o} {input.a} {input.r} {input.t} {input.p} {output.flag} true;
		"""

rule score_motifs:
	input:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		p = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_PROMOTERS"]]),
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],""]),
		b = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"motif_inference/network_bins/"]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.infer_motifs"])
	output:
		m = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"motif_inference/motifs.txt"]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.score_motifs"])
	shell:
		"""
		printf "\e[1m[Step 5.1 completed]\e[0m\n\n\e[1mStep 5.2: Scoring promoters with inferred motifs\e[0m\n"; ./CODE/run_score_motifs.sh {input.o} {input.b} {input.r} {input.p} {output.m} {output.flag} true;
		"""

rule build_motif_network:
	input:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		g = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_GENES"]]),
		i = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"motif_inference/motifs.txt"]),
		f = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"motif_inference/motifs_score/"]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.score_motifs"])
	params:
		v = config["MOTIF_THRESHOLD"]
	output:
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/mn.adjmtr"])
	shell:
		"""
		printf "\e[1m[Step 5.2 completed]\e[0m\n\n\e[1mStep 5.3: Mapping motif network\e[0m\n"; python CODE/build_motif_network.py -i {input.i} -r {input.r} -g {input.g} -f {input.f} -t robust -v {params.v} -o {output.o}; printf "\e[1m[Step 5.3 completed]\e[0m\n\n";
		"""

rule assemble_final_network:
	input:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		a = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["DBD_PID_DIR"]]),
		d = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],"networks/"]),
		i = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/npwa_bnwa.adjmtr"]),
		m = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/mn.adjmtr"])
	output:
		o = NETPROPHET2_NETWORK
	shell:
		"""
		printf "\e[1mStep 6: Assemble final network\e[0m\n"; python CODE/combine_networks.py -s resort -n {input.i} -b {input.m} -od {input.d} -om npwa_bnwa_mn.adjmtr; python CODE/weighted_avg_similar_dbds.py -n {input.d}/npwa_bnwa_mn.adjmtr -r {input.r} -a {input.a} -d 50 -f single_dbds -o {output.o}; printf "\e[1m[Step 6 completed]\e[0m\n\n\e[1mNetProphet 2.0 network is ready\e[0m at %s\n" {output.o};
		"""

