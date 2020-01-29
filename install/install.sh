#!/usr/bin/env bash
# Install all dependencies

#########################################################
## System dependencies

# tools
SYSTEM_DEPS="Rscript R make cmake python3 python"

# libs
SYSTEM_PACKS="pandoc"

## Software to install locally
## TOOLS
ALL_SOFT="blast R_packages"


## Versions and URLs
blast_VERSION=2.9.0
blast_URL=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-${blast_VERSION}+-x64-linux.tar.gz

function install_blast {
	pinfo "Installing blast..."
	pushd $TEMP_FOLDER
	wget -c $blast_URL -O tmp.tar.gz
	tar zxvpf tmp.tar.gz
	cp ncbi-blast-${blast_VERSION}+/bin/* $INSTALL_BIN
	## taxonomy
	wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
	mkdir -p $INSTALL_DIR/db
	mv taxdb.tar.gz $INSTALL_DIR/db
	pushd $INSTALL_DIR/db
	tar xzvf taxdb.tar.gz
	rm -f taxdb.tar.gz
	popd
	popd
	pinfo "Installing blast...done."
}


function install_R_packages {
    mkdir -p $INSTALL_DIR/Rlibs
    R_LIBS_USER=$INSTALL_DIR/Rlibs R --vanilla <<EOF
repo<-"http://www.stats.bris.ac.uk/R/"

########################
# Check if version is ok
version <- getRversion()
currentVersion <- sprintf("%d.%d", version\$major, version\$minor)
message("R version:",version)
usebiocmanager<-TRUE
if ( version\$major < 3 || (version\$major>=3 && version\$minor<5) ) {
  cat("ERROR: R version should be 3.5 or above\n")
  q(status=1)
}

########################
# Where to install the packages
assign(".lib.loc",.libPaths()[1],envir=environment(.libPaths))

message("Using library: ", .libPaths()[1])
##print(.libPaths())

message("_____________________________________________________")

if (version\$major > 3 || (version\$major == 3 && version\$minor>5)) {
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager",repo=repo)
   BiocManager::install()
} else {
   usebiocmanager<-FALSE
   source("http://bioconductor.org/biocLite.R")
} 

message("_____________________________________________________")

message("Installing packages")
## R PACKAGES to install
packages2install<-c("data.table","devtools")

for (p in packages2install ) {
  message("PACKAGE:",p,"\n")
  if ( usebiocmanager ) BiocManager::install(p,ask=FALSE)
  else  biocLite(p,ask=FALSE)
}

#message("PACKAGE:","d3treeR","\n")
#devtools::install_github("timelyportfolio/d3treeR")

EOF

    
}


######################################################################
function pinfo {
    echo "[INFO] $*"
}

function usage {
    echo "Usage: install.sh  [-x all|tool_name -i toplevel_installation_folder]  ";
    echo " -x software: install/update software.";
    echo " -i dir : install/update all files to directory 'dir' (default: $PWD/bastools)";
}

function check_system_deps {
    local bin
    pinfo "Checking dependencies..."
    local MISSING=0
    for bin in $SYSTEM_DEPS; do
	local PATH2BIN=`which $bin 2> /dev/null`
	if [ "$PATH2BIN-" == "-" ]; then
	    pinfo " $bin not found!"
	    #
	    MISSING=1
	else
	    pinfo " $bin found: $PATH2BIN"
	fi
    done
    pinfo "Checking dependencies...done."
    if [ $MISSING == 1 ]; then
	pinfo "ERROR: Unable to proceed"
	exit 1
    fi

}

function install_all {
    check_system_deps    
    for tt in $ALL_TOOLS; do
	install_$tt
    done
}

## default installation folder
INSTALL_DIR=$PWD/bastools
set +eux
if [ "$BASTOOLS_DIR-" != "-" ]; then
    ## update previous installation
    INSTALL_DIR=$BASTOOLS_DIR
fi

## by default install all software
MODE=all
DEBUG=0
set +u
while getopts "i:x:hH"  Option
do
    case $Option in
	i) INSTALL_DIR=$OPTARG;;
	x) MODE=$OPTARG;;
	d) DEBUG=1;;
	h ) usage; exit;;
	H ) usage; exit;;
    esac
done

if [ $DEBUG == 1 ]; then
    set -eux
else
    set -eu
fi


BASTOOLS_DIR=$INSTALL_DIR

INSTALL_DIR=$(readlink -f $INSTALL_DIR)
## PREFIX=INSTALL_DIR
INSTALL_BIN=$INSTALL_DIR/bin
TEMP_FOLDER=$INSTALL_DIR/tmp
set -eu

mkdir -p $INSTALL_BIN
TEMP_FOLDER=$(mktemp -d)
##-p $PWD)
mkdir -p $TEMP_FOLDER

function gen_env_sh {
    env_file=$1
    pinfo "Generating $1..."
    python_dir=python$(python --version 2> /dev/stdout | sed "s/.* \(.*\)\..*/\1/")
    python3_dir=python$(python3 --version 2> /dev/stdout | sed "s/.* \(.*\)\..*/\1/")
    set +u
    if [ "$PYTHONPATH-" == "-" ]; then
	PYTHONPATH=.
    fi
    set -u
    cat <<EOF > $env_file
# source $env_file
# Created `date`
export PATH=$INSTALL_BIN:\$PATH    
export BASTOOLS_DIR=$INSTALL_DIR
set +u
## Python
export PYTHONPATH=$BASTOOLS_DIR/lib64/$python_dir/site-packages:\$BASTOOLS_DIR/lib/$python_dir/site-packages:\$BASTOOLS_DIR/lib64/$python3_dir/site-packages:\$BASTOOLS_DIR/lib/$python3_dir/site-packages:$PYTHONPATH
## R packages
export R_LIBS=\$BASTOOLS_DIR/Rlibs:\$R_LIBS

# BLAST
if [ "\$BLASTDB-" == "-" ]; then
export BLASTDB=\$BASTOOLS_DIR/db
fi
EOF
    pinfo "Generating $1...done"
}

BASTOOLS_ENV_FILE=$INSTALL_DIR/msi_env.sh
if [ ! -e $BASTOOLS_ENV_FILE ]; then
    gen_env_sh $BASTOOLS_ENV_FILE    
fi

source $BASTOOLS_ENV_FILE

x="-$(echo $ALL_SOFT make|sed -E 's/\s+/-|-/g')|-all-"
echo $x
if [[ "-$MODE-" =~ ^($x) ]]; then
    pinfo "Installation mode: $MODE"
else
    pinfo Valid values for -x: $ALL_TOOLS
    echo ERROR: invalid value $MODE for -x parameter
    rm -rf $TEMP_FOLDER
    exit 1
fi

set -u
pinfo "Installation folder: $INSTALL_DIR"
install_$MODE
rm -rf $TEMP_FOLDER
exit 0
