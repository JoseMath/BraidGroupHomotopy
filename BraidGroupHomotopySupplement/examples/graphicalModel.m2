

A=matrix{
    {0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1},
    {0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1},
    {0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1},
    {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1},
    {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1},
    {0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1},
    {0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1},
    {0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1}
    }
xList= for i to numrows  A-1 list value("x"|i)
kk=ZZ/30103
R=kk[xList|{s}]**kk[t]
xList= for i to numrows  A-1 list value("x"|i)
#xList
F=sum for a in entries transpose A list product apply(a,xList,(i,j)->j^i)
uList=for i to numcols A-1 list random(1,30103)
uList=(for i to numcols A-1 list random(1,30103))+t*for i to numcols A-1 list random(1,30103)
lhs=flatten entries (A*transpose matrix{uList})
rhs= uList//sum*s*for i in xList list (i*diff(i,F))

eqs={1-s*F}|lhs-rhs
I=ideal eqs
det matrix makeJac(eqs,{s}|xList)


