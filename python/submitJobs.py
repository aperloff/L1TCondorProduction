from jobSubmitterL1T import jobSubmitterL1T

def submitJobs():  
    mySubmitter = jobSubmitterL1T()
    mySubmitter.run()
    
if __name__=="__main__":
    submitJobs()
