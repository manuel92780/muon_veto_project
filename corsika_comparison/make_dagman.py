
#arguement parser
import argparse
parser = argparse.ArgumentParser(description='Make dagman for corsika processing')
parser.add_argument('--did', default=20694, type=int, help='Simpord Number for this sim')
args = parser.parse_args()
    
f = open('dagman_'+str(args.did)+'.dag','w')
file_name = 'corsika_'+str(args.did)+'.txt'
file_list = open(file_name).read().splitlines()
n_files = len(file_list)

print "will submit " + str(n_files) + " jobs"
for njob in range(n_files):
    a = "JOB job"+str(njob)+" job.sub\n"
    b = "VARS job"+str(njob)+" njobs=\""+str(n_files) +"\" did=\""+str(args.did)+"\" seed=\""+str(njob)+"\" file_name=\""+file_list[njob]+"\"\n"
    f.write(a)
    f.write(b)

f.close()
