
f = open('dagman.dag','w')

file_name = 'corsika_list.txt'
file_list = open(file_name).read().splitlines()
n_files = 1
n_files = len(file_list)

print "will submit " + str(n_files) + " jobs"
for njob in range(n_files):
    did = 0
    if '20209' in file_list[njob]:
        did = 20209
    if '20243' in file_list[njob]:
        did = 20243  
    a = "JOB job"+str(njob)+" job.sub\n"
    b = "VARS job"+str(njob)+" did=\""+str(did)+"\" seed=\""+str(njob)+"\" file_name=\""+file_list[njob]+"\"\n"
    f.write(a)
    f.write(b)

f.close()
