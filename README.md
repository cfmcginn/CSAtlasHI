# CSAtlasHI

On Lxplus or ACF/RCF, recommend the following setup commands:

setupATLAS\
lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt pythia8 240" #PYTHIA8\
lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt fastjet" #FASTJET\
lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt fjcontrib" #FASTJET CONTRIB\
lsetup "root 6.18.04-x86_64-centos7-gcc8-opt" #root6 distribution

on a personal machine you will have to build these packages and link manually

If running a macro complains 'error while loading shared libraries: libCSATLAS.so', try from the directory of the Makefile:\
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib
