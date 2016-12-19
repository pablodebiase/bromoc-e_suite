#!/usr/bin/python
import sys

def istitle(string):
    head=['ATOM','BOND','ANGL','THET','DIHE','IMPR','IMPH','NONB','NBON','CMAP','NBFI','HBON','PHI','END']
    for word in head:
        if string.startswith(word):
            return True
    return False

def rmnum(string):
    st=''
    for wrd in string.split():
        try:
            a=float(wrd)
        except:
            st=st+' '+wrd
    return st.strip()

def main():
    if len(sys.argv) < 3:
        print 'This program reduces the parameters to the minimum needed'
        print 'USAGE: ',sys.argv[0], ' file.psf file.prm'
        return
    
    with open(sys.argv[1], 'r') as psf:
        # Read first
        line=psf.readline().strip().split()
        if line[0] != 'PSF':
            print 'Not a PSF'
            return
        psf.readline()
        # Read title
        line=psf.readline().strip().split()
        for i in range(int(line[0])):
             psf.readline().strip()
        psf.readline()
        # Read NAtom
        line=psf.readline().strip().split()
        spec=set()
        for i in range(int(line[0])):
             spec.add(psf.readline().strip().split()[5])
        psf.close()
    if len(spec) == 0:
        return
    spec.add('X')
    virtual='' 
    with open(sys.argv[2], 'r') as prm:
        for line in prm.readlines():
            pl=line.split('!')[0].strip()
            if len(pl) == 0:
                #prmout.write(line)  # uncomment to allow single comment lines 
                continue
            if len(pl.split('*')[0].strip()) == 0:
                #prmout.write(line) 
                virtual+=line
                continue 
            if istitle(pl):
                #prmout.write(line)
                virtual+=line
                continue
            labels=rmnum(pl).split()
            if len(labels) > 0 and set(labels).issubset(spec):
                #prmout.write(line)
                virtual+=line
                continue
    outname='out.prm'
    prmout=open(outname, 'w')
    pt=False
    pline=None
    for line in virtual.splitlines():
        if (istitle(line.split('!')[0].strip())):
            if pline != None and not pt:
                prmout.write(pline+'\n')
            pt=True
        else:
            if pline != None:
                prmout.write(pline+'\n')
            pt=False
        pline=line
    prmout.close()
    print outname, 'was saved'
    return

if __name__ == "__main__":
    main()

