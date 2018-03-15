
restart
A=matrix{
    {0,0,0,0,1,1,1,1},
    {0,0,1,1,0,0,1,1},
    {0,1,0,1,0,1,0,1},
    {0,0,0,0,0,0,1,1},
    {0,0,0,1,0,0,0,1},
    {0,0,0,0,0,1,0,1}}

xList=for i from 1 to numrows A list value("x"|i)
uList=for i from 1 to 8 list value("u"|i)
kk=QQ
R=kk[{s}|xList]**kk[uList]

xList=for i from 1 to numrows A list value("x"|i)
uList=for i from 1 to 8 list value("u"|i)
mons=apply(entries transpose A,i->product apply(i,xList,(a,b)->b^a))

F=sum mons
eq1={1-s*F}
eq2=flatten entries((sub(A,R)*transpose matrix{uList})-transpose matrix{apply(xList,i->((sum uList)*s*i*diff(i,F)))})
EQ=eq1|eq2

degree sub(ideal EQ,apply(uList,u->u=>random(1,30103)))
E=eliminate({s}|drop(xList,1),ideal EQ)
--discriminant(E_0,x1)



