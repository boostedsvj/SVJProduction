# SVJProduction

## GEN-SIM production

To make GEN-SIM samples, `CMSSW_7_1_26` is used.
All of the necessary setup (including installation of the latest Pythia8 version, CMSSW package updates and compilation, and checkout of this repo)
is performed by [setup.sh](./scripts/setup.sh):
```
wget https://raw.githubusercontent.com/kpedro88/SVJProduction/master/scripts/setup.sh
chmod +x setup.sh
./setup.sh -c CMSSW_7_1_26 -p
cd CMSSW_7_1_26/src
cmsenv
```

To run the sample production interactively with example parameters:
```
cd SVJ/Production/test
cmsRun runSVJ.py config=SVJ.Production.step1_GEN output=RAWSIMoutput outpre=step1 mZprime=2000.0 mDark=20.0 rinv=0.3 alpha=0.1 part=1 maxEvents=10
```

## cmsDriver commands

GEN only:
```
cmsDriver.py SVJ.Production.EmptyFragment_cff --python_filename step1_GEN.py --mc --eventcontent RAWSIM --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --datatier GEN-SIM --conditions MCRUN2_71_V3::All --beamspot Realistic50ns13TeVCollision --step GEN --magField 38T_PostLS1 --fileout file:step1.root --no_exec
```