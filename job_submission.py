import submission
import shutil
import shelve
import os
import subprocess
import nbody

### SELECT ONLY THE BODIES YOU WANT TO SUBMIT ###
body_list = ['one_body','two_body','three_body','four_body','five_body','six_body','seven_body']
##body_list += ['two_body','three_body']
#body_list = ['four_body']
#body_list = ['five_body','six_body']

### PASS n_body_target (MAX NUMBER OF n-BODIES FOR THIS CALCULATION) ##
n_body_target = 7
nbody.check_child_job_status(n_body_target)

#for body in body_list:
#    status = shelve.open('{0}/status'.format(body))
#    for job in status['job_status']:
#        dir = submission.get_dir(job)
#        os.chdir('{0}/{1}'.format(body,dir))
#        shutil.copy('../../submit','submit_{0}'.format(dir))
#        submit_dir = 'submit_{0}'.format(dir)
#        subprocess.call(['qsub', submit_dir])
#        os.chdir('../..')
#    status.close()

#for body in body_list:
#    status = shelve.open('{0}/status'.format(body))
#    job_list = status['job_status']
#    print job_list
#    while job_list: 
#        job,val = job_list.popitem(last=True)
#        dir = submission.get_dir(job)
#        os.chdir('{0}/{1}'.format(body,dir))
#        shutil.copy('../../submit','submit_{0}'.format(dir))
##        print dir 
#        submit_dir = 'submit_{0}'.format(dir)
#        subprocess.call(['qsub', submit_dir])
#        os.chdir('../..')
#    status.close()

for body in body_list:
    status = shelve.open('{0}/status'.format(body))
    job_list = status['job_status']
    print job_list
'''    for key in job_list: 
	if job_list[key] == "not_started" or job_list[key] == "error": 
#	if job_list[key] == "not_started": # Use this line of you want to skip over "error" calculations for a bit 
            #job,val = job_list.popitem(last=True)
            dir = submission.get_dir(key)
            os.chdir('{0}/{1}'.format(body,dir))
            shutil.copy('../../submit','submit_{0}'.format(dir))
            #print dir 
            submit_dir = 'submit_{0}'.format(dir)
            subprocess.call(['qsub', submit_dir])
            os.chdir('../..')
    status.close()
'''       
nbody.check_child_job_status(n_body_target)
