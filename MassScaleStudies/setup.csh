echo "rootsys: "$ROOTSYS
echo "setenv ROOTSYS /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.00/x86_64-slc5-gcc43-opt/root/"
setenv ROOTSYS /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.00/x86_64-slc5-gcc43-opt/root/
echo "rootsys: "$ROOTSYS
setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:$LD_LIBRARY_PATH
setenv PATH ${ROOTSYS}/bin:$PATH

setenv XRDCP /afs/cern.ch/sw/lcg/external/xrootd/3.1.0p2/x86_64-slc5-gcc43-opt/
source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.csh

setenv LD_LIBRARY_PATH `pwd`/lib:$XRDCP/lib64:$LD_LIBRARY_PATH
setenv PATH `pwd`/bin/:$XRDCP/bin/:$PATH
