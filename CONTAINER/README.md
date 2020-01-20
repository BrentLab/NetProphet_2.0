# NETPROPHET 2.0 CONTAINER

We use [Singularity](https://sylabs.io/singularity/) to run NetProphet in a container environment.  
Provided is the definition file that allows you to easily build a distributable NetProphet image.

Some useful applications of running in a container:

    • Allows you to run customized software applications in your own environment  
    • Allows you to run an application on a cluster without actually installing anything 
    • Allows you to run a series of applications (a 'pipeline') built on different platforms or environments
    • Distribute your pipeline/workflow to others 

This guide will show you how to build a NetProphet container using Singularity and Conda.


### SYSTEM REQUIREMENTS

* Singularity (>= v3.5.2, tested on v3.5.2)
* Cloned copy of NetProphet_2.0 repository


### BUILDING INSTRUCTIONS

1. Install Singularity
	
    Follow SyLabs instructions to install Singularity:
    [Singularity Installation](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps)

2. Build the Image:

    You can choose to build an image locally (more CPU intensive):

	```
    cd NetProphet_2.0/CONTAINER
    singularity build --sandbox np2image.sif np2.def
	```

    Or remotely using Sylabs container services:
    Visit [Sylabs Cloud](https://cloud.sylabs.io/auth) to create an access token, then:

	```
    cd NetProphet_2.0/CONTAINER
    singularity build --remote --sandbox np2image.sif np2.def
	```
    
    This will produce a SIF image file. 


### RUNNING CONTAINER

1. Example usage:

	```
    singularity run -B <path to NetProphet_2.0>:/NetProphet_2.0 np2image.sif NetProphet2 -s -f config.json
	```

    This will attach the local NetProphet_2.0 directory to the image at location /NetProphet_2.0 and execute NetProphet2 in serial mode.

    


