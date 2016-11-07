#!/usr/bin/python
# This program creates effective potentials file from the output of bromoc efpot write function
# The output files are compatible to start Inverse Monte Carlo iteration process from imc-macro
from os import path
def addtolist(mylist,element):
    if not element in mylist:
        mylist.append(element)
    return 
strt=0.1
ends=25.0
step=0.1
types=[]
finp=open('efpot.dat','r')
fout=None
for line in finp:
    if ' - ' in line:
        if fout:
            fout.close() 
        
        name=line.translate(None,' \n\t\r').lower()
        addtolist(types,name.split('-')[0])
        addtolist(types,name.split('-')[1])
        fname=name+'.pot'
        fout=open(fname,'w')
    else:
        x, y = map(float,line.split())[0:2]
        num=x/step
        inum=int(num)
        diff=abs(inum-num)
        if diff < 1e-8 and x > strt-step and x < ends+step:
            fout.write('%16.10f     %25.18e\n' % (inum*step, y))
finp.close()
if fout:
    fout.close()
types=sorted(types)
ntyp=len(types)
fout=open('system-in.pot','w')
fout.write('%12i%12i%13.8f     %13.8f\n' % (ntyp, int(ends/step), 0.0, ends))
for i in range(ntyp):
    for j in range(i,ntyp):
        name=types[i]+'-'+types[j]+'.pot'
        exist=False
        if path.isfile(name):
            exist=True
        else:
            name=types[j]+'-'+types[i]+'.pot'
            if path.isfile(name):
                exist=True
        if exist:
            finp=open(name,'r')
            for line in finp:
                fout.write('%s%12i%12i\n' % (line.translate(None,'\n\t\r'), i+1, j+1))
            finp.close()
fout.close()
