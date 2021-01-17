
p_src_code = "/scratch/mblab/dabid/netprophet/code_netprophet2.1/"
p_lib = "/home/dabid/R/x86_64-pc-linux-gnu-library/3.4"

# install lars library, tested for R 3.4.3
install.packages(paste(p_src_code, "SRC/R_pkgs/lars_1.2.tar.gz", sep="")
                , repos=NULL
                , type="source"
                , lib=p_lib
                )
 
# isntall Rmpi ( the open mpi have to be installed or loaded)
R CMD INSTALL /scratch/mblab/dabid/netprophet/code_netprophet2.1/SRC/R_pkgs/Rmpi_0.6-3.tar.gz --configure-args="--with-Rmpi-include=/opt/apps/openmpi/1.8.8/include  --with-Rmpi-type=OPENMPI --with-Rmpi-libpath=/opt/apps/openmpi/1.8.8/lib/"

# install Bayes trees
install.packages("/scratch/mblab/dabid/netprophet/code_netprophet2.1/SRC/R_pkgs/BayesTree_0.3-1.3.tar.gz", repos=NULL, type="source")

# install abind
install.packages("abind")
# then select mirror
