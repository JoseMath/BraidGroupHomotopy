

restart
loadPackage"BraidGroupHomotopy"
--installPackage"Bertini"
matLabTriangle2=(nB,aTri)->(
    print ("x"|nB|"Loop = "|toString new Array from  ((SetDownstairsStartPoint|aTri|SetDownstairsStartPoint)/realPart));
print ("y"|nB|"Loop= "|toString new Array from  ((SetDownstairsStartPoint|aTri|SetDownstairsStartPoint)/imaginaryPart)
    ))
matLabTL2=(n,nB,aTri)->(
    print ("t"|nB|"x"|n|" = "|toString new Array from  (drop(aTri,-1)/realPart));
print ("t"|nB|"y"|n|"= "|toString new Array from  (drop(aTri,-1)/imaginaryPart)
    ))

--
printingPrecision=300
R=CC[z,t,x,y,yA,s]
varList=(z,t,x,y,yA,s)
beta=(34-6*sqrt(21))/7
f=z^3*(z^2-2*z+beta)^2-t/20
computeBranchPoints(z,t,f)--there should be two branch points. 
#SetBraidGroupBranchPoints===2
SetDownstairsStartPoint={.05-0.05*ii}
SetUpstairsStartFiber=bertiniZeroDimSolve({sub(f,{t=>first SetDownstairsStartPoint})},AffVariableGroup=>{z}    )/coordinates//flatten//sort
#SetDownstairsStartPoint==1
#SetUpstairsStartFiber==first degree f
SetEncirclingTriangles=computeEncirclingTriangles(null)
--Displays loops for matlab
apply(#SetEncirclingTriangles,i->matLabTriangle2(i+1,SetEncirclingTriangles_i))
theBraids=apply(SetEncirclingTriangles,oneTriangle->computeBraid(varList,f,oneTriangle))


-----Matlab output for first branch point. 
matLabTriangle2(1,first SetEncirclingTriangles)
tw1 = {.0569997066067544-.0495679396844149*ii, .0920432315249522-.0474048609604063*ii}
tw2 = {.102323604262132-.0383226120528186*ii, .145701708096762-3.02535774210355e-15*ii, .156945590491976+.00993346653975299*ii}
tw3 = {.154830219478869+.0106429369753943*ii, .144915950278069+.0139680656858538*ii, .132069655589161+.018276561090605*ii, .124453856678343+.020830810032281*ii, .0748379627828145+.0374713944206531*ii}
tw4 = {.0768283471651006+.0276525200008866*ii, .0768753291922017+.0274207503848279*ii, .0824337904763846-9.0205620750794e-17*ii, .0920432315249522-.0474048609604063*ii}
tw5 = {.0851842569159786-.0478282345245617*ii, .0535297429062691-.0497821248918425*ii, .0509466401696213-.0499415681722949*ii, .05-.05*ii}
twList1={tw1,tw2,tw3,tw4,tw5}
apply(#twList1,i->matLabTL2(i+1,1,twList1_i))

--Matlab output for second branchPoint
tw1 = {.0498445352738245-.0406513393541629*ii, .0496848414004356-.0310483646713677*ii, .0491688532319607-.0000201153508513768*ii, .0491685187212098-3.46944695195361e-17*ii, .0490033288919836+.00993346653976056*ii}
tw2 = {.0414107829811976+.0124799166998389*ii, .040977618999594+.0126251947795839*ii, .00995567810080619+.0230295871001387*ii, -.0331042988171783+.0374713944206607*ii}
tw3 = {-.0308597588615528+.026398731195117*ii, -.0289196049344247+.0168276513881125*ii, -.0255803179919943+.000354431729081506*ii, -.0255084711236076-1.12410081243297e-15*ii, -.0255082584545641-.00000104912932106377*ii, -.0248115370145241-.00343808389818381*ii, -.0230146454648722-.0123024282087949*ii, -.0223759721019042-.0154531028216831*ii, -.0210060880460832-.0222109530331743*ii, -.0158990300750407-.0474048609603987*ii}
tw4 = {-.00932015958754227-.041592723518908*ii, -.00845543784004356-.0408287805814905*ii, .0105296051976244-.0240563453014192*ii, .0377567348781195-.00000239559363188896*ii, .0377594464967627+1.17961196366423e-16*ii, .0490033288919836+.00993346653976056*ii}
tw5 = {.0490716475908192+.00582521416064698*ii, .0491123751812949+.0033761156953807*ii, .0491685187212096+5.20417042793042e-18*ii, .049347228388074-.0107464636560645*ii, .0499319861024346-.0459100764322451*ii, .05-.05*ii}
twList2={tw1,tw2,tw3,tw4,tw5}
matLabTriangle2(2,last SetEncirclingTriangles)
apply(#twList2,i->matLabTL2(i+1,2,twList2_i))



{*--Output

(Twist locus on a segment, (.05-.05*ii, .0920432315249522-.0474048609604063*ii), {.0569997066067544-.0495679396844149*ii, .0920432315249522-.0474048609604063*ii})
tw1 = {.0569997066067544-.0495679396844149*ii, .0920432315249522-.0474048609604063*ii}
(Twist locus on a segment, (.0920432315249522-.0474048609604063*ii, .156945590491976+.00993346653975299*ii), {.102323604262132-.0383226120528186*ii, .145701708096762-3.02535774210355e-15*ii, .156945590491976+.00993346653975299*ii})
tw2 = {.102323604262132-.0383226120528186*ii, .145701708096762-3.02535774210355e-15*ii, .156945590491976+.00993346653975299*ii}
(Twist locus on a segment, (.156945590491976+.00993346653975299*ii, .0748379627828145+.0374713944206531*ii), {.154830219478869+.0106429369753943*ii, .144915950278069+.0139680656858538*ii, .132069655589161+.018276561090605*ii, .124453856678343+.020830810032281*ii, .0748379627828145+.0374713944206531*ii})
tw3 = {.154830219478869+.0106429369753943*ii, .144915950278069+.0139680656858538*ii, .132069655589161+.018276561090605*ii, .124453856678343+.020830810032281*ii, .0748379627828145+.0374713944206531*ii}
(Twist locus on a segment, (.0748379627828145+.0374713944206531*ii, .0920432315249522-.0474048609604063*ii), {.0768283471651006+.0276525200008866*ii, .0768753291922017+.0274207503848279*ii, .0824337904763846-9.0205620750794e-17*ii, .0920432315249522-.0474048609604063*ii})
tw4 = {.0768283471651006+.0276525200008866*ii, .0768753291922017+.0274207503848279*ii, .0824337904763846-9.0205620750794e-17*ii, .0920432315249522-.0474048609604063*ii}
(Twist locus on a segment, (.0920432315249522-.0474048609604063*ii, .05-.05*ii), {.0851842569159786-.0478282345245617*ii, .0535297429062691-.0497821248918425*ii, .0509466401696213-.0499415681722949*ii, .05-.05*ii})
tw5 = {.0851842569159786-.0478282345245617*ii, .0535297429062691-.0497821248918425*ii, .0509466401696213-.0499415681722949*ii, .05-.05*ii}
Twist loci found! 
(END SEGMENT 0, .0569997066067544-.0495679396844149*ii)
(END SEGMENT 1, .0920432315249522-.0474048609604063*ii)
(END SEGMENT 2, .102323604262132-.0383226120528186*ii)
The twist is : [1,2]
The twist is : [3,4]
The twist is : [5,6]
3 twists on fiber over: 
.1457
{-.11066-.13341*ii, -.11066+.13341*ii, .37864-.1269*ii, .37864+.1269*ii, 1.0553-.09867*ii, 1.0553+.09867*ii, 1.3535}
{-.11066+.13341*ii, -.11066-.13341*ii, .37864+.1269*ii, .37864-.1269*ii, 1.0553+.09867*ii, 1.0553-.09867*ii, 1.3535}
(END SEGMENT 3, .145701708096762-3.02535774210355e-15*ii)
(END SEGMENT 4, .156945590491976+.00993346653975299*ii)
(END SEGMENT 5, .154830219478869+.0106429369753943*ii)
(END SEGMENT 6, .144915950278069+.0139680656858538*ii)
(END SEGMENT 7, .132069655589161+.018276561090605*ii)
(END SEGMENT 8, .124453856678343+.020830810032281*ii)
(END SEGMENT 9, .0748379627828145+.0374713944206531*ii)
(END SEGMENT 10, .0768283471651006+.0276525200008866*ii)
(END SEGMENT 11, .0768753291922017+.0274207503848279*ii)
The twist is : [2,1]
1 twists on fiber over: 
.082434
{-.091732+.11521*ii, -.091732-.11521*ii, .26599, .49478, .95519, 1.1323, 1.3352}
{-.091732-.11521*ii, -.091732+.11521*ii, .26599, .49478, .95519, 1.1323, 1.3352}
(END SEGMENT 12, .0824337904763846-9.0205620750794e-17*ii)
(END SEGMENT 13, .0920432315249522-.0474048609604063*ii)
(END SEGMENT 14, .0851842569159786-.0478282345245617*ii)
(END SEGMENT 15, .0535297429062691-.0497821248918425*ii)
(END SEGMENT 16, .0509466401696213-.0499415681722949*ii)
(END SEGMENT 17, .05-.05*ii)
END BRANCH POINT 

(Twist locus on a segment, (.05-.05*ii, .0490033288919836+.00993346653976056*ii), {.0498445352738245-.0406513393541629*ii, .0496848414004356-.0310483646713677*ii, .0491688532319607-.0000201153508513768*ii, .0491685187212098-3.46944695195361e-17*ii, .0490033288919836+.00993346653976056*ii})
tw1 = {.0498445352738245-.0406513393541629*ii, .0496848414004356-.0310483646713677*ii, .0491688532319607-.0000201153508513768*ii, .0491685187212098-3.46944695195361e-17*ii, .0490033288919836+.00993346653976056*ii}
(Twist locus on a segment, (.0490033288919836+.00993346653976056*ii, -.0331042988171783+.0374713944206607*ii), {.0414107829811976+.0124799166998389*ii, .040977618999594+.0126251947795839*ii, .00995567810080619+.0230295871001387*ii, -.0331042988171783+.0374713944206607*ii})
tw2 = {.0414107829811976+.0124799166998389*ii, .040977618999594+.0126251947795839*ii, .00995567810080619+.0230295871001387*ii, -.0331042988171783+.0374713944206607*ii}
(Twist locus on a segment, (-.0331042988171783+.0374713944206607*ii, -.0158990300750407-.0474048609603987*ii), {-.0308597588615528+.026398731195117*ii, -.0289196049344247+.0168276513881125*ii, -.0255803179919943+.000354431729081506*ii, -.0255084711236076-1.12410081243297e-15*ii, -.0255082584545641-.00000104912932106377*ii, -.0248115370145241-.00343808389818381*ii, -.0230146454648722-.0123024282087949*ii, -.0223759721019042-.0154531028216831*ii, -.0210060880460832-.0222109530331743*ii, -.0158990300750407-.0474048609603987*ii})
tw3 = {-.0308597588615528+.026398731195117*ii, -.0289196049344247+.0168276513881125*ii, -.0255803179919943+.000354431729081506*ii, -.0255084711236076-1.12410081243297e-15*ii, -.0255082584545641-.00000104912932106377*ii, -.0248115370145241-.00343808389818381*ii, -.0230146454648722-.0123024282087949*ii, -.0223759721019042-.0154531028216831*ii, -.0210060880460832-.0222109530331743*ii, -.0158990300750407-.0474048609603987*ii}
(Twist locus on a segment, (-.0158990300750407-.0474048609603987*ii, .0490033288919836+.00993346653976056*ii), {-.00932015958754227-.041592723518908*ii, -.00845543784004356-.0408287805814905*ii, .0105296051976244-.0240563453014192*ii, .0377567348781195-.00000239559363188896*ii, .0377594464967627+1.17961196366423e-16*ii, .0490033288919836+.00993346653976056*ii})
tw4 = {-.00932015958754227-.041592723518908*ii, -.00845543784004356-.0408287805814905*ii, .0105296051976244-.0240563453014192*ii, .0377567348781195-.00000239559363188896*ii, .0377594464967627+1.17961196366423e-16*ii, .0490033288919836+.00993346653976056*ii}
(Twist locus on a segment, (.0490033288919836+.00993346653976056*ii, .05-.05*ii), {.0490716475908192+.00582521416064698*ii, .0491123751812949+.0033761156953807*ii, .0491685187212096+5.20417042793042e-18*ii, .049347228388074-.0107464636560645*ii, .0499319861024346-.0459100764322451*ii, .05-.05*ii})
tw5 = {.0490716475908192+.00582521416064698*ii, .0491123751812949+.0033761156953807*ii, .0491685187212096+5.20417042793042e-18*ii, .049347228388074-.0107464636560645*ii, .0499319861024346-.0459100764322451*ii, .05-.05*ii}
Twist loci found! 
(END SEGMENT 0, .0498445352738245-.0406513393541629*ii)
(END SEGMENT 1, .0496848414004356-.0310483646713677*ii)
(END SEGMENT 2, .0491688532319607-.0000201153508513768*ii)
The twist is : [1,2]
1 twists on fiber over: 
.049169
{-.077153-.10036*ii, -.077153+.10036*ii, .19449, .56803, .89118, 1.179, 1.3216}
{-.077153+.10036*ii, -.077153-.10036*ii, .19449, .56803, .89118, 1.179, 1.3216}
(END SEGMENT 3, .0491685187212098-3.46944695195361e-17*ii)
(END SEGMENT 4, .0490033288919836+.00993346653976056*ii)
(END SEGMENT 5, .0414107829811976+.0124799166998389*ii)
(END SEGMENT 6, .040977618999594+.0126251947795839*ii)
(END SEGMENT 7, .00995567810080619+.0230295871001387*ii)
(END SEGMENT 8, -.0331042988171783+.0374713944206607*ii)
(END SEGMENT 9, -.0308597588615528+.026398731195117*ii)
(END SEGMENT 10, -.0289196049344247+.0168276513881125*ii)
(END SEGMENT 11, -.0255803179919943+.000354431729081506*ii)
The twist is : [2,3]
The twist is : [4,5]
The twist is : [6,7]
3 twists on fiber over: 
-.025508
{-.099487, .041668-.11222*ii, .041668+.11222*ii, .7357-.10321*ii, .7357+.10321*ii, 1.2724-.04557*ii, 1.2724+.04557*ii}
{-.099487, .041668+.11222*ii, .041668-.11222*ii, .7357+.10321*ii, .7357-.10321*ii, 1.2724+.04557*ii, 1.2724-.04557*ii}
(END SEGMENT 12, -.0255084711236076-1.12410081243297e-15*ii)
(END SEGMENT 13, -.0255082584545641-.00000104912932106377*ii)
(END SEGMENT 14, -.0248115370145241-.00343808389818381*ii)
(END SEGMENT 15, -.0230146454648722-.0123024282087949*ii)
(END SEGMENT 16, -.0223759721019042-.0154531028216831*ii)
(END SEGMENT 17, -.0210060880460832-.0222109530331743*ii)
(END SEGMENT 18, -.0158990300750407-.0474048609603987*ii)
(END SEGMENT 19, -.00932015958754227-.041592723518908*ii)
(END SEGMENT 20, -.00845543784004356-.0408287805814905*ii)
(END SEGMENT 21, .0105296051976244-.0240563453014192*ii)
(END SEGMENT 22, .0377567348781195-.00000239559363188896*ii)
The twist is : [1,2]
1 twists on fiber over: 
.037759
{-.070562-.093376*ii, -.070562+.093376*ii, .17047, .59251, .86905, 1.1935, 1.3156}
{-.070562+.093376*ii, -.070562-.093376*ii, .17047, .59251, .86905, 1.1935, 1.3156}
(END SEGMENT 23, .0377594464967627+1.17961196366423e-16*ii)
(END SEGMENT 24, .0490033288919836+.00993346653976056*ii)
(END SEGMENT 25, .0490716475908192+.00582521416064698*ii)
(END SEGMENT 26, .0491123751812949+.0033761156953807*ii)
The twist is : [2,1]
1 twists on fiber over: 
.049169
{-.077153+.10036*ii, -.077153-.10036*ii, .19449, .56803, .89118, 1.179, 1.3216}
{-.077153-.10036*ii, -.077153+.10036*ii, .19449, .56803, .89118, 1.179, 1.3216}
(END SEGMENT 27, .0491685187212096+5.20417042793042e-18*ii)
(END SEGMENT 28, .049347228388074-.0107464636560645*ii)
(END SEGMENT 29, .0499319861024346-.0459100764322451*ii)
(END SEGMENT 30, .05-.05*ii)
END BRANCH POINT 

a_1-a_3-a_5 a_1^{-1}
a_1 a_2-a_4-a_6 a_1 a_1^{-1}


*}