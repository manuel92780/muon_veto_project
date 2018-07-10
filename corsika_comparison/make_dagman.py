
f = open('dagman.dag','w')

file_name = 'corsika_list.txt'
file_list = open(file_name).read().splitlines()
#n_files = len(file_list)
n_files = 1

print "will submit " + str(n_files) + " jobs"
for njob in range(n_files):
    a = "JOB job"+str(njob)+" job.sub\n"
    b = "VARS job"+str(njob)+" seed=\""+str(njob)+"\" file_name=\""+file_list[njob]+"\"\n"
    f.write(a)
    f.write(b)

f.close()
