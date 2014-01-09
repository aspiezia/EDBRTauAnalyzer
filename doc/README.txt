####### TAKE THE CODE FORM GIT #######
cmsrel CMSSW_5_3_13
cd CMSSW_5_3_13
cmsenv
git cms-merge-topic -u cms-tau-pog:CMSSW_5_3_X_boostedTaus
cd src
cmsenv
git clone git@github.com:YOURNAME/EDBRTauAnalyzer.git
setenv CVSROOT ":ext:YOURNAME@lxplus5.cern.ch:/afs/cern.ch/user/c/cvscmssw/public/CMSSW" 
cvs co -r V00-02-03s TauAnalysis/CandidateTools
scram b


###### RUN ON THE QUEUES #####
#1) Inside data/ create a txt file with name dataset_fileList.txt and put the list of root file that you want to process
#2) In createConfigFile.sh correct the python config file that you want to use
#3) Run on the queues with the following command:
./runOnQueue.sh dataset N queue TypeOfAnalysis NameOfFolder Site DATA_or_MC SR_or_sideband

#where:
# - dataset            is the name of the process your are studying that shoud match data/dataset_fileList.txt
# - N                  gives the number of jobs to run through the following formula N(jobs) = N(root file)/N
# - queue              is the queue on which you want to run (unused if you are running at Zurich)
# - TypeOfAnalysis     can be "fullyLeptonic" or "semiLeptonic"
# - NameOfFolder       is the name of the forlder in which you are saving the output
# - Site               can be CERN or ZURICH depending where you are running
# - data_or_MC         tells if you are processing data or MC
# - SR_or_sideband:    SR for the Signal Region and sideband for the sideband in the jet mass