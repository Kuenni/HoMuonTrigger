#!/usr/bin/env bash

#Add the python directory to the path and the Pythonpath
BASEDIR="$( cd "$( dirname "$BASH_SOURCE[0]" )" && pwd )"
echo "$BASEDIR"
export HOMUONTRIGGER_BASE=$BASEDIR
export PATH=$HOMUONTRIGGER_BASE/python:$PATH
export PYTHONPATH=$HOMUONTRIGGER_BASE/python:$PYTHONPATH

submitData(){
	runParallel.py --conditions=data -n 339  --sourceFiles files_SingleMuon_Run2015D-PromptReco-v4_RECO --lumiFile Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt
}

collectData(){
	runParallel.py -c -d /user/kuensken/tapasTasks/collisionData2015 -s 10 -o SingleMuon2015D.root
}

submitNoPu(){
	runParallel.py --conditions=noPu -n 200
}

collectNoPu(){
	runParallel.py -c -d /user/kuensken/tapasTasks/SingleMuonGunPooja -s 10 -o NoPUAnalyzed.root
}

submitPu(){
	runParallel.py --conditions=pu -n 400 
}

collectPu(){
	runParallel.py -c -d /user/kuensken/tapasTasks/SingleMuonGunWithPU -s 10 -o SingleMuWithPu.root
}

export -f submitData
export -f collectData
export -f submitNoPu
export -f collectNoPu
export -f submitPu
export -f collectPu
