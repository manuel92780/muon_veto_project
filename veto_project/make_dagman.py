
f = open('dagman.dag','w')

step   = 10000
nevents = 1000001

print "will submit " + str((nevents-1)/step) + " jobs"
for njob in range(step, nevents, step):
    start = njob-step+1
    end   = njob
    a = "JOB job"+str(njob/step)+" job.sub\n"
    b = "VARS job"+str(njob/step)+" nevents=\""+str(step)+"\" start=\""+str(start)+"\" end=\""+str(end)+"\"\n"
    f.write(a)
    f.write(b)

f.close()
