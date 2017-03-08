configfile: "config.json"
workdir: config["NETPROPHET2_DIR"]

NETPROPHET2_NETWORK = config["FILENAME_NETPROPHET2_NETWORK"]

rule all:
	input:
		NETPROPHET2_NETWORK

rule make_directories:
	output:
		r = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp",
		n = config["NETPROPHET2_DIR"] + "/OUTPUT/networks",
		m = config["NETPROPHET2_DIR"] + "/OUTPUT/motif_inference",
		s = config["NETPROPHET2_DIR"] + "/OUTPUT/motif_inference/network_scores/",
		b = config["NETPROPHET2_DIR"] + "/OUTPUT/motif_inference/network_bins/",
		p = config["NETPROPHET2_DIR"] + "/OUTPUT/motif_inference/motifs_pfm/",
		q = config["NETPROPHET2_DIR"] + "/OUTPUT/motif_inference/motifs_score/",
		flag = config["NETPROPHET2_DIR"] + "/LOG/flag.make_dir"
	shell:
		"""
		mkdir -p {output.r}; mkdir -p {output.n}; mkdir -p {output.m}; mkdir -p {output.s}; \
		mkdir -p {output.b}; mkdir -p {output.p}; mkdir -p {output.q}; touch {output.flag}
		"""

rule prepare_resources:
	input:
		g = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_GENES"],
		r = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_REGULATORS"],
		e = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_EXPRESSION_DATA"],
		f = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_FOLDCHANGE_DATA"],
		c = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_SAMPLE_CONDITIONS"],
		flag = config["NETPROPHET2_DIR"] + "/LOG/flag.make_dir"
	output:
		r = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/rdata.expr",
		f = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/data.fc.tsv",
		a = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/allowed.adj",
		p1 = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/data.pert.adj",
		p2 = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/data.pert.tsv",
		l = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/regulator_sublists/",
		flag = config["NETPROPHET2_DIR"] + "/LOG/flag.prepare_resources"
	shell:
		"""
		python CODE/prepare_resources.py -g {input.g} -r {input.r} -e {input.e} \
		-f {input.f} -c {input.c} -or {output.r} -of {output.f} -oa {output.a} \
		-op1 {output.p1} -op2 {output.p2} -ol {output.l}; touch {output.flag}
		"""

rule map_np_network:
	input:
		g = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_GENES"],
		f = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_REGULATORS"],
		t = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_EXPRESSION_DATA"],
		d = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_DE_ADJMTR"],
		u = config["NETPROPHET2_DIR"] + "/SRC/NetProphet1/",
		r = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/rdata.expr",
		a = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/allowed.adj",
		p = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/data.pert.adj",
		flag = config["NETPROPHET2_DIR"] + "/LOG/flag.prepare_resources"
	output:
		o = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/",
		n = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/np.adjmtr"
	shell:
		"""
		./SRC/NetProphet1/netprophet -m -c -u {input.u} -t {input.t} -r {input.r} \
		-a {input.a} -p {input.p} -d {input.d} -g {input.g} -f {input.f} \
		-o {output.o} -n {output.n}
		"""

rule map_bart_network:
	input:
		f = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_REGULATORS"],
		t = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/data.fc.tsv",
		p = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/data.pert.tsv",
		flag = config["NETPROPHET2_DIR"] + "/LOG/flag.prepare_resources"
	output:
		o = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/bn.adjmtr"
	shell:
		"""
		sbatch CODE/run_build_bart_network.sh {input.t} {input.p} {input.f} {output.o}
		"""

rule weighted_average_np_network:
	input:
		r = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_REGULATORS"],
		a = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["DIR_DBD_PID"],
		n = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/np.adjmtr"
	output:
		o = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/npwa.adjmtr"
	shell:
		"""
		python CODE/weighted_avg_similar_dbds.py -n {input.n} -r {input.r} \
		-a {input.a} \-d 50 -t single_dbds -o {output.o}
		"""

rule weighted_average_bart_network:
	input:
		r = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_REGULATORS"],
		a = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["DIR_DBD_PID"],
		n = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/bn.adjmtr"
	output:
		o = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/bnwa.adjmtr"
	shell:
		"""
		python CODE/weighted_avg_similar_dbds.py -n {input.n} -r {input.r} \
		-a {input.a} -d 50 -t single_dbds -o {output.o}
		"""

rule combine_npwa_bnwa:
	input:
		n = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/npwa.adjmtr",
		b = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/bnwa.adjmtr"
	output:
		o = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/npwa_bnwa.adjmtr"
	shell:
		"""
		Rscript CODE/quantile_combine_networks.r {input.n} {input.b} {output.o}
		"""

rule infer_motifs:
	input:
		r = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_REGULATORS"],
		t = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_GENES"],
		p = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_PROMOTERS"],
		o = config["NETPROPHET2_DIR"] + "/OUTPUT/",
		a = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/npwa_bnwa.adjmtr",
		l = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/regulator_sublists/"
	params:
		m = config["PROMOTER_LENGTH"]
	output:
		flag = config["NETPROPHET2_DIR"] + "/LOG/flag.infer_motifs"
	shell:
		"""
		bash CODE/run_infer_motifs.sh {input.o} {input.a} {input.r} {input.t} \
		{input.l} {input.p} {params.m} {output.flag}
		"""

rule score_motifs:
	input:
		r = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_REGULATORS"],
		p = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_PROMOTERS"],
		o = config["NETPROPHET2_DIR"] + "/OUTPUT/",
		b = config["NETPROPHET2_DIR"] + "/OUTPUT/motif_inference/network_bins/",
		l = config["NETPROPHET2_DIR"] + "/RESOURCES/tmp/regulator_sublists/",
		flag = config["NETPROPHET2_DIR"] + "/LOG/flag.infer_motifs"
	output:
		m = config["NETPROPHET2_DIR"] + "/OUTPUT/motif_inference/motifs.txt",
		flag = config["NETPROPHET2_DIR"] + "/LOG/flag.score_motifs"
	shell:
		"""bash CODE/run_score_motifs.sh {input.o} {input.b} {input.r} {input.l} \
		{input.p} {output.m} {output.flag}
		"""

rule build_motif_network:
	input:
		r = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_REGULATORS"],
		g = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_GENES"],
		i = config["NETPROPHET2_DIR"] + "/OUTPUT/motif_inference/motifs.txt",
		f = config["NETPROPHET2_DIR"] + "/OUTPUT/motif_inference/motifs_score/",
		flag = config["NETPROPHET2_DIR"] + "/LOG/flag.score_motifs"
	params:
		v = config["MOTIF_THRESHOLD"]
	output:
		o = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/mn.adjmtr"
	shell:
		"""
		python CODE/build_motif_network.py -i {input.i} -r {input.r} -g {input.g} \
		-f {input.f} -t robust -v {params.v} -o {output.o}
		"""

rule assemble_final_network:
	input:
		r = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["FILENAME_REGULATORS"],
		a = config["NETPROPHET2_DIR"] + "/RESOURCES/" + config["DIR_DBD_PID"],
		d = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/",
		i = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/npwa_bnwa.adjmtr",
		m = config["NETPROPHET2_DIR"] + "/OUTPUT/networks/mn.adjmtr"
	output:
		o = NETPROPHET2_NETWORK
	shell:
		"""python CODE/combine_networks.py -s resort -n {input.i} -b {input.m} \
		-od {input.d} -om npwa_bnwa_mn.adjmtr; \
		python CODE/weighted_avg_similar_dbds.py -n {input.d}/npwa_bnwa_mn.adjmtr \
		-r {input.r} -a {input.a} -d 50 -f single_dbds -o {output.o}; \
		rm LOG/*; \
		rm -rf OUTPUT/motif_inference; \
		rm -rf OUTPUT/networks; \
		rm -rf RESOURCES/tmp; \
		echo '### COMPLETED! ###';
		"""

