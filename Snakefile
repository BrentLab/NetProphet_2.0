configfile: "config.json"
workdir: config["NETPROPHET2_DIR"]

NETPROPHET2_NETWORK = "/".join([config["OUTPUT_DIR"],
							config["FILENAME_NETPROPHET2_NETWORK"]])

rule all:
	input:
		NETPROPHET2_NETWORK

rule make_directories:
	output:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],"tmp"]),
		n = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],"networks"]),
		m = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],"motif_inference"]),
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
		mkdir -p {output.r}; mkdir -p {output.n}; mkdir -p {output.m}; mkdir -p {output.s}; \
		mkdir -p {output.b}; mkdir -p {output.p}; mkdir -p {output.q}; touch {output.flag}
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
		l = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/regulator_sublists/"]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.prepare_resources"])
	shell:
		"""
		python CODE/prepare_resources.py -g {input.g} -r {input.r} -e {input.e} \
		-c {input.c} -or {output.r} -of {output.f} -oa {output.a} \
		-op1 {output.p1} -op2 {output.p2} -ol {output.l}; touch {output.flag}
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
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.prepare_resources"])
	output:
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/"]),
		n = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"np.adjmtr"])
	shell:
		"""
		./SRC/NetProphet1/netprophet -m -c -u {input.u} -t {input.t} -r {input.r} \
		-a {input.a} -p {input.p} -d {input.d} -g {input.g} -f {input.f} \
		-o {output.o} -n {output.n}
		"""

rule map_bart_network:
	input:
		f = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		t = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/data.fc.tsv"]),
		p = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/data.pert.tsv"]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.prepare_resources"])
	output:
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/bn.adjmtr"])
	shell:
		"""
		sbatch CODE/run_build_bart_network.sh {input.t} {input.p} {input.f} {output.o}
		"""

rule weighted_average_np_network:
	input:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		a = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["DIR_DBD_PID"]]),
		n = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/np.adjmtr"])
	output:
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/npwa.adjmtr"])
	shell:
		"""
		python CODE/weighted_avg_similar_dbds.py -n {input.n} -r {input.r} \
		-a {input.a} \-d 50 -t single_dbds -o {output.o}
		"""

rule weighted_average_bart_network:
	input:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		a = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["DIR_DBD_PID"]]),
		n = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/bn.adjmtr"])
	output:
		o = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/bnwa.adjmtr"])
	shell:
		"""
		python CODE/weighted_avg_similar_dbds.py -n {input.n} -r {input.r} \
		-a {input.a} -d 50 -t single_dbds -o {output.o}
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
		Rscript CODE/quantile_combine_networks.r {input.n} {input.b} {output.o}
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
					"networks/npwa_bnwa.adjmtr"]),
		l = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/regulator_sublists/"])
	params:
		m = config["PROMOTER_LENGTH"]
	output:
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.infer_motifs"])
	shell:
		"""
		bash CODE/run_infer_motifs.sh {input.o} {input.a} {input.r} {input.t} \
		{input.l} {input.p} {params.m} {output.flag}
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
		l = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					"tmp/regulator_sublists/"]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.infer_motifs"])
	output:
		m = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"motif_inference/motifs.txt"]),
		flag = "/".join([config["NETPROPHET2_DIR"],"LOG/flag.score_motifs"])
	shell:
		"""
		bash CODE/run_score_motifs.sh {input.o} {input.b} {input.r} {input.l} \
		{input.p} {output.m} {output.flag}
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
		python CODE/build_motif_network.py -i {input.i} -r {input.r} -g {input.g} \
		-f {input.f} -t robust -v {params.v} -o {output.o}
		"""

rule assemble_final_network:
	input:
		r = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["FILENAME_REGULATORS"]]),
		a = "/".join([config["NETPROPHET2_DIR"],config["RESOURCES_DIR"],
					config["DIR_DBD_PID"]]),
		d = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],"networks/"]),
		i = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/npwa_bnwa.adjmtr"]),
		m = "/".join([config["NETPROPHET2_DIR"],config["OUTPUT_DIR"],
					"networks/mn.adjmtr"]),
		f = config["CLEANUP"]
	output:
		o = NETPROPHET2_NETWORK
	shell:
		"""
		python CODE/combine_networks.py -s resort -n {input.i} -b {input.m} \
		-od {input.d} -om npwa_bnwa_mn.adjmtr; \
		python CODE/weighted_avg_similar_dbds.py -n {input.d}/npwa_bnwa_mn.adjmtr \
		-r {input.r} -a {input.a} -d 50 -f single_dbds -o {output.o}; \
		rm LOG/flag.*
		if {input.cleanup}; then
			rm LOG/*; 
			rm -rf OUTPUT/motif_inference; 
			rm -rf OUTPUT/networks; 
			rm -rf RESOURCES/tmp;
		fi
		echo '### COMPLETED! ###';
		"""

