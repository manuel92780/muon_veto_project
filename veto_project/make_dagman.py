
f = open('dagman.dag','w')

step   = 100000
nevents = 8000001
lowE=1e3
midE=1e4
highE=1e7
njobs = (nevents-1)/step

#do lowE first
print "will submit " + str(njobs) + " jobs"
for njob in range(step, nevents, step):
    start = njob-step+1
    end   = njob
    a = "JOB job"+str(njob/step)+" job.sub\n"
    b = "VARS job"+str(njob/step)+" njobs=\""+str(njobs) +"\" lowE=\""+str(lowE)+"\" highE=\""+str(highE)+"\" nevents=\""+str(step)+"\" start=\""+str(start)+"\" end=\""+str(end)+"\"\n"
    f.write(a)
    f.write(b)

#do highE second
#for njob in range(nevents, 2*nevents-1, step):
#    start = njob
#    end   = njob+step
#    a = "JOB job"+str(njob/step + 1)+" job.sub\n"
#    b = "VARS job"+str(njob/step +1)+" lowE=\""+str(midE)+"\" highE=\""+str(highE)+"\" nevents=\""+str(step)+"\" start=\""+str(start)+"\" end=\""+str(end)+"\"\n"
#    f.write(a)
#    f.write(b)


f.close()
