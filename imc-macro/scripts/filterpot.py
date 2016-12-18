#!/usr/bin/python
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob,copy
from scipy.fftpack import fft,ifft,fftshift,ifftshift

def genplot(xdat,y1dat,y2dat=[],y3dat=[],y4dat=[],y5dat=[],title=''):
    x=xdat
    y1=y1dat
    plt.plot(x,y1,'k')
    if len(y2dat) > 1:
        y2=y2dat
        plt.plot(x,y2,'r')
    if len(y3dat) > 1:
        y3=y3dat
        plt.plot(x,y3,'g')
    if len(y4dat) > 1:
        y4=y4dat
        plt.plot(x,y4,'b')
    if len(y5dat) > 1:
        y5=y5dat
        plt.plot(x,y5,'y')
    plt.title(title)
    plt.grid()
    plt.show()
    return

def readfile(filename):
    x=[]
    y=[]
    with open(filename,'r') as xy:
        for ln in xy.readlines():
            xi, yi = map(float,ln.split()[0:2])
            x.append(xi)
            y.append(yi)
    return np.array(x), np.array(y)

def sync(L, dec):
    l=int(L/2)
    F=0.5/dec
    b=2.0*F*np.sinc(2.0*F*np.linspace(-l,l,L))
    bsum=np.sum(b)
    b/=bsum
    #genplot(np.linspace(-l,l,L),b)
    return b

def fsync(data,dec):
    n=len(data)
    # find minimum L for optimum result
    #a=1.0 # bad
    #a=2  # bad
    #a=2.45902475 #not too bad
    #a=3.47088974 #not too bad
    a=1.43029666 #good truncation
    L=int(2.0*dec*a)
    if L % 2 == 0:
        L+=1
    b = sync(L, dec)
    i=0
    lm=int(L)/2
    cdata=data.copy()
    while i<n:
        cdata[i]=data[i]*b[lm]
        l=1
        while l <= lm:
            if i+l < n:
                cdata[i]+=data[i+l]*b[lm+l]
            else:
                cdata[i]+=data[2*n-i-l-2]*b[lm+l]
            if i-l >= 0:
                cdata[i]+=data[i-l]*b[lm-l]
            else:
                cdata[i]+=(2.0*data[i]-data[l-i])*b[lm-l]
            l+=1
        i+=1
    return cdata

def fwindow(data,dec):
    n=len(data)
    #Size of extra head/tail
    n25=int(n*0.5)
    #Decimated size
    m=(n+2*n25)/dec/2
    #Add extra head
    tdata=np.fliplr([data[1:n25+1]])[0]
    idata=np.append(tdata,data)
    #Add extra tail
    tdata=np.fliplr([data[-n25-1:n-1]])[0]
    idata=np.append(idata,tdata)
    #plot3d(np.linspace(-0.125,1.125,n+2*n25), idata)
    #FFT
    fdata=fft(idata)
    #Kill high frequencies
    fdata[m:-m]=0
    #IFFT
    odata=ifft(fdata)
    #Chop off extra tail and head
    odata=odata[n25:n+n25]
    #genplot(np.linspace(-1,1,n), np.abs(fftshift(fft(data))),np.abs(fftshift(fft(odata))))
    return odata

def writepot(fname, x, y):
    out=open(fname, 'w')
    for i in range(len(x)):
        out.write('%16.8f      %16.8f\n' % (x[i], y[i]))
    out.close()
    return

# Single potential file
def main2():
    if len(sys.argv) < 2:
        print 'No arguments given'
        return
    filename=sys.argv[1]
    if len(sys.argv) < 3:
        dec=3.0
    else:
        dec=float(sys.argv[2])
    x,y = readfile(filename)
    yf=fsync(y, dec)
    writepot('out.pot', x, yf)
    genplot(x,y,yf)
    return
    km=1.0/(x[1]-x[0])
    k=np.linspace(-km*0.5,km*0.5,len(y))
    kk=np.linspace(-len(y)*0.5,len(y)*0.5,len(y))
    iy=fftshift(fft(y))
    genplot(k,iy)
    return

def getdict(filename):
    rdfdict={}
    with open(filename,'r') as rdf:
        spec=int(rdf.readline().split()[0])
        for i in range(spec):
            rdf.readline()
        for line in rdf.readlines():
            xi, yi = map(float,line.split()[0:2])
            it, jt = map(int,line.split()[2:4])
            if yi == 0.0:
               rdfdict['%d,%d' % (min(it,jt),max(it,jt))] = xi
    return rdfdict

            
# Multiple potential file
def main():
    if len(sys.argv) < 2:
        print 'No arguments given'
        return
    filename=sys.argv[1]
    if len(sys.argv) < 3:
        dec=3.0
    else:
        dec=float(sys.argv[2])
    inp=open(filename,'r')
    out=open('out.pot','w')
    rdfs=glob.glob('*.rdf')
    isrdf=len(rdfs) > 0
    rdf=''
    if isrdf:
        rdfdict=getdict(rdfs[0])
    print rdfdict
    line=inp.readline()
    out.write(line)
    species, points = map(int,line.split()[0:2])
    potentials=species*(species+1)/2
    for i in range(potentials):
        x=[]
        y=[]
        for j in range(points):
            line=inp.readline()
            xi, yi = map(float,line.split()[0:2])
            it, jt = map(int,line.split()[2:4])
            if isrdf:
                fac=rdfdict['%d,%d' % (min(it,jt),max(it,jt))]
            else:
                fac=3.0
            if isrdf and xi <= fac:
                out.write('%21.16f %23.17G %12i %12i\n' % (xi,yi,it,jt))
                
            elif not isrdf and yi > 21.0:
                out.write('%21.16f %23.17G %12i %12i\n' % (xi,yi,it,jt))
            else:
                x.append(xi)
                y.append(yi)
        x=np.array(x)
        y=np.array(y)
        i0=0 #np.abs(x-fac-1.0).argmin()+1
        i1=np.abs(x-5.0).argmin()+1
        i2=np.abs(x-10.0).argmin()+1
        i3=np.abs(x-15.0).argmin()+1
        yf=copy.copy(y)
        yf[i0:i1+1]=fsync(y, 3.0)[i0:i1+1]
        yf[i1:i2+1]=fsync(y, 5.0)[i1:i2+1]
        yf[i2:i3+1]=fsync(y, 7.0)[i2:i3+1]
        yf[i3:]=fsync(y, 9.0)[i3:]
        for j in range(len(x)):
            out.write('%21.16f %23.17G %12i %12i     %23.17G\n' % (x[j],yf[j],it,jt,y[j]))
        print it,jt
        genplot(x,y,yf)
    inp.close()
    out.close()
    return

if __name__ == "__main__":
    if 'multi' in sys.argv[0]:
        main()
    else:
        main2()

