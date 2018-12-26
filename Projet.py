from numpy import *
from math import *
from matplotlib.pyplot import *
from numpy.linalg import *
from cmath import *
def argsinh(x):
    return(log(x+sqrt(x**2+1)))
def f(x):
    return(2*(x**2)*(1-cosh((B[0]-A[0])/x))+l**2-(A[1]-B[1])**2)
def df(x):
    return(4*x*(1-cosh((B[0]-A[0])/x))+2*(B[0]-A[0])*sinh((B[0]-A[0])/x))
def verif(A,B,l):
    if ((A[0]==B[0])or(l<=(sqrt((A[0]-B[0])**2+(A[1]-B[1])**2)))):
        return(0)
    else:
        return(1)
#********************************************tache1********************************************
def tache1(A,B,l):
    if (verif(A,B,l)==0):
        print("Attention verifier que les conditions la chainette n'existe pas!!!")
        return(0,0,0,0)
    else:
        print("la chainette existe!!!!!")
        x=0.4
        y=x-(f(x)/df(x))
        N_iter=0
        while ((abs(y-x)>10**(-5))and(N_iter<10000)and(df(x)!=0)):
            N_iter+=1
            x=y
            y=x-(f(x)/df(x))
        alpha=y
        beta=((((A[0]+B[0])/alpha)-argsinh((A[1]-B[1])/(2*alpha*sinh((A[0]-B[0])/2))))*0.5*alpha)
        gam=A[1]-alpha*cosh((A[0]-beta)/alpha)
        return(alpha,beta,gam,N_iter)
def chainette(alpha,beta,gam,x):
    return(alpha*cosh((x-beta)/alpha)+gam)
def tracer(alpha,beta,gam):
    X=[]
    Y=[]
    h=abs(B[0]-A[0])/100
    j=min(A[0],B[0])
    n=0;
    while(n<=100):
        n+=1
        X.append(j)
        Y.append(chainette(alpha,beta,gam,j))
        j+=h
    plot([A[0],A[0]],[0,j-h])
    plot([B[0],B[0]],[0,chainette(alpha,beta,gam,j-h)])
    plot(X,Y)
    title("Traçage de la chainette ")
    text(A[0],A[1],'A')
    text(B[0],B[1],'B')
    show()
#********************************************tache2********************************************
#********************************************tache3********************************************
def Jh(Y):
   s1=0
   s2=0
   x1=xa
   x2=x1+h
   for i in range(len(Y)-1):
      y=Y[i]+((((x2-x1)/h)*(Y[i+1]-Y[i])))
      s1+=y*sqrt(1+((Y[i+1]-Y[i])/h)**2)
      s2+=sqrt(1+((Y[i+1]-Y[i])/h)**2)
      x1=x2
      x2=x2+h
   return(h*s1+(1/eps)*((h*s2-l)**2))
#********************************************tache4********************************************
def goldensearch(a,b,f):
    eps=10**(-5)
    erreur=abs(b-a)
    t=(1+sqrt(5))/2
    iterations=0
    while((erreur)>eps)and(iterations<nbritermax):
        a1=real(a+((b-a)/(t**2)))
        b1=real(b-((b-a)/(t**2)))
        if (f(a1)>f(b1)):
            a=a1
        elif (f(a1)<f(b1)):
            b=b1
        elif (f(a1)==f(b1)):
            a=a1
            b=b1
        erreur=b-a
        iterations=iterations+1      
    return(a,iterations)
#********************************************tache5********************************************
def comp_grad(j,Y):
    Y=[real(k) for k in Y]
    N=size(Y)
    s=0
    for i in range(1,N-1):
        s+=(1+((Y[i+1]-Y[i])/h)**2)**0.5
    s=s-l
    s=s*(2*(h**2))/eps
    c3=((Y[j+1]-Y[j])/h)**2
    c4=((Y[j]-Y[j-1])/h)**2
    c1=(1+c4)**0.5
    c2=(1+c3)**0.5
    s=s*((c4/c1)-(c3/c2))
    c5=(X[j+1]-X[j])/2
    c6=(X[j]-X[j-1])/2
    s=s+h-c5*(c2+c3/c2)+c6*(c1+c4/c1)
    return(s)
#********************************************tache6********************************************
def DJh(Y):
    N=size(Y)
    grad=zeros(N)
    for j in range(1,N-1):
        grad[j]=comp_grad(j,Y)
    return( grad )
#********************************************tache7********************************************
def res(Y):
    while(True):
        f = lambda r: Jh(Y-(dot(r,DJh(Y))))
        r=goldensearch(0,999,f)[0]
        X=Y
        Y=Y-(dot(r,DJh(Y)))
        if(norm(X-Y)<10**(-4)):
            break
    return(Y)
#********************************************tache8********************************************
#Ceci n'est pas complet
#********************************************tache9********************************************
#********************************************programme principal********************************************
print("donner les points d'attache et la longeur")
A=[0,0]
B=[0,0]
l=0
A[0]=float(input('Donner xa:'))
A[1]=float(input('Donner ya:'))
B[0]=float(input('Donner xb:'))
B[1]=float(input('Donner yb:'))
l=float(input('Donner la longeur l:'))
alpha,beta,gam,n=tache1(A,B,l)
print("alpha=",alpha,"beta=",beta,"gamma=",gam,"nbre itération",n)
tracer(alpha,beta,gam)

nbritermax=100
xa=0
ya=2
xb=2
yb=3
l=10
eps=10**(-5)
N=1000
h=(xb-xa)/(N-1)

f=lambda x : cosh(x-(1/2))-cosh(1/2)
N=int((xb-xa)/h)
X=[i*h for i in range(N)]
Y=[f(i)+2 for i in X]
Y[0]=ya
Y[N-1]=yb
Y=res(Y)
Y[0]=ya
Y[N-1]=yb
plt.plot(X,Y)
plt.show()
