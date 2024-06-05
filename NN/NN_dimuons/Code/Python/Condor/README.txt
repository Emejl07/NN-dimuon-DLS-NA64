Training the NN requires some workload that should not be run on /eos/

This is why Condor is used. A very simple implementation of Condor has been done here for simply running the DNN_dimuons_pyTorch.py file in one job. 

As of today, one should not place logs mentioned in the HTCondor submit file to EOS. This means that this repository should be in another folder than in /eos/. For my case, /afs/cern.ch/work was used.

To initialize the run, select the end destination of the model and the plots in DNN_dimuons_pyTorch.py and also your environment in run_job.sh and condor.submit aswell as the location for the log,output and error file. 

After this is done, which by the way should be mostly in /afs/, read more about this here to avoid any error: https://batchdocs.web.cern.ch/troubleshooting/eos.html#no-eos-submission-allowed,
simply run with the command: condor_submit condor.submit

Check run with: condor_q
Check output during runtime with: condor_tail

  
