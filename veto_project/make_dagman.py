
f = open('dagman.dag','w')

step   = 10000
nevents = 1000001
lowE=1e1
midE=1e4
highE=1e6

#do lowE first
print "will submit " + str(2*(nevents-1)/step) + " jobs"
for njob in range(step, nevents, step):
    start = njob-step+1
    end   = njob
    a = "JOB job"+str(njob/step)+" job.sub\n"
    b = "VARS job"+str(njob/step)+" lowE=\""+str(lowE)+"\" highE=\""+str(midE)+"\" nevents=\""+str(step)+"\" start=\""+str(start)+"\" end=\""+str(end)+"\"\n"
    f.write(a)
    f.write(b)

#do highE second
for njob in range(nevents, 2*nevents-1, step):
    start = njob
    end   = njob+step
    a = "JOB job"+str(njob/step + 1)+" job.sub\n"
    b = "VARS job"+str(njob/step +1)+" lowE=\""+str(midE)+"\" highE=\""+str(highE)+"\" nevents=\""+str(step)+"\" start=\""+str(start)+"\" end=\""+str(end)+"\"\n"
    f.write(a)
    f.write(b)

f.close()
