# TrackingPerf

//----------------------Paul Vaucelle----------------------//
//     Ph.D student for CMS, @IPHC, 10/2022->10/2025       //
//              DisplacedTop analysis code                 //
// contact me @: paul.vaucelle@iphc.cnrs.fr or             //
//               paul.vaucelle@cern.ch                     //
//---------------------------------------------------------//

The analysis code is available and can run on RECO data tier.
A MiniAOD version is available @: https://github.com/Threshic/FlyingTop/tree/master/plugins
For the RECO Version : 
Look for the _RECO version or "Version-05_10_2022-AVF (Pour présentation TPS)" of TrackingPerf.cc (in plugins) on the history version of the repo.
Finally, since some data is lost due to MINIAOD format, a comparison is made for the firsthit in the tracker between RECO and MINIAOD on the version 
of 16/11/2022.
+ Don't pay attention to the _MINIAOD since these files are a mix of RECO and MINIAOD have a high chanc eof crashing when running on MINIAOD files

For people @IPHC that want access to the code:

# Release 10_6_20_FLY : 

###
### création du répertoire :
###

cd /opt/sbg/cms/<insert uiX_dataY>/<insert name>

export V0_CMS_SW_DIR=/cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/crab3/crab.sh 
export SCRAM_ARCH=slc7_amd64_gcc700

scramv1 p CMSSW CMSSW_10_6_20
mv CMSSW_10_6_20 CMSSW_10_6_20_LLP
cd CMSSW_10_6_20_LLP
scramv1 b ProjectRename
cd src
eval  `scramv1 r -sh`

git cms-init

# cp -r /opt/sbg/cms/ui2_data1/blochd/CMSSW_10_6_20_FLY_pourPaul/src/TrackingPerf . 
cp -r /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/TrackingPerf .

# just for information for the very first time:
# mkdir TrackingPerf
# cd TrackingPerf
# mkedanlzr TrackingPerf
# x TrackingPerf/plugins/BuildFile.xml & // and add <use name="DataFormats/PatCandidates"/>

###
### compilation :
###

scramv1 b clean
scramv1 b -j6 

-----------------------------------

###
### exécution :
###

cd /opt/sbg/cms/<insert uiX_dataY>/<insert name>/CMSSW_10_6_20_LLP/src

export V0_CMS_SW_DIR=/cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700

eval `scramv1 runtime -sh`

export PYTHONPATH=$PYTHONPATH:/opt/glite/lib64/python:/opt/fpconst/lib/python2.4/site-packages
export PYTHONPATH=${PYTHONPATH}:${GLITE_LOCATION}/lib64
cmsenv

# pour specifier une version de CRAB:
source /cvmfs/cms.cern.ch/crab3/crab.sh 
voms-proxy-init -rfc -voms cms -valid 192:00

cd FlyingTop/FlyingTop/test
cmsRun step3_2018.py
# crab submit -c crab_config_mc_2018.py
