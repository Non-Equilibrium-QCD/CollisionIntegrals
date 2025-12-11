p1,p2,q,w=var('p1 p2 q w')
assume(p1>0,p2>0,q>0)

s=solve([p1>1/2*(q+w),p2>1/2*(q-w),-q<w, w<q,p1>0,p2>0,q>0],w)
print(s)

s=solve([p1>1/2*(q+w),p2>1/2*(q-w),-q<w, w<q,p1>0,p2>0,q>0],q)
print(s)


c1,s1,cq,phi1q=var('c1 s1 cq phi1q')

solve([s1^2 * (1-cq*cq) * sin(phi1q)^2 == (x-c1*cq)^2], cq)
