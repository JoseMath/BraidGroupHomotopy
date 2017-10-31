
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
alpha=(34+6*sqrt(21))/7
f=-1/1000*z^3*(z^2-2*z+alpha)^2-t

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


matLabTriangle2(1,first SetEncirclingTriangles)

-----Matlab output for first branch point. 
tw1 = {.139278220184272-.0476175796315843*ii, .147249586032148-.047404860960407*ii}
tw2 = {.153497010396114-.0418855413252677*ii, .157739843322658-.0381371884505414*ii, .175200289648387-.0227116660846958*ii, .182992922881545-.0158272236157688*ii, .200561373056886-.000306284688375305*ii, .200660054056451-.000219104444680547*ii, .20066448058917-.000215193801265938*ii, .200717611456599-.000168255059354951*ii, .200908062603961+9.57567358739198e-16*ii, .212151944999172+.0099334665397523*ii}
tw3 = {.203579515259291+.0128085581545913*ii, .199192077034142+.0142800530789608*ii, .198863231264089+.0143903440635172*ii, .192433236038032+.0165468884644459*ii, .13004431729001+.0374713944206524*ii}
tw4 = {.133964753976512+.018131273066208*ii, .13764014498358-1.80411241501588e-16*ii, .140548235762022-.0143460622023175*ii, .144292991207221-.0328195207847156*ii, .146420755919154-.0433161136474567*ii, .147249586032148-.047404860960407*ii}
tw5 = {.139372777271377-.0476150563427887*ii, .0620967876846337-.0496771930117645*ii, .05-.05*ii}
twList1={tw1,tw2,tw3,tw4,tw5}
apply(#twList1,i->matLabTL2(i+1,1,twList1_i))

--Matlab output for second branchPoint
tw1 = {.0499649775475425-.0478939722785278*ii, .0495940983990646-.0255916572433402*ii, .0491878892647111-.00116482138012693*ii, .0491685187212788+1.52655665885959e-16*ii, .0490136282555203+.00931412827396011*ii, .0490033288920667+.00993346653975239*ii}
tw2 = {.0297825501079181+.0163798885682567*ii, .0214277296070334+.019181996616556*ii, .0169808654119673+.0206734223086301*ii, .00515923208561742+.0246382584033228*ii, .00283970560009898+.02541620017125*ii, .00201954931854355+.0256912708940327*ii, -.00631591506102611+.0284868871269537*ii, -.0119704556044553+.0303833531939715*ii, -.0160530497270709+.0317526070092857*ii, -.0256515708655207+.0349718375515459*ii, -.0263694747059806+.0352126140145329*ii, -.0331042988170953+.0374713944206525*ii}
tw3 = {-.0306720609536044+.0254727883457755*ii, -.0306588803499729+.0254077663867164*ii, -.0286613257154425+.0155535200021312*ii, -.0285600422432051+.015053872947246*ii, -.0255084711235156-1.76941794549634e-15*ii, -.0240126386326793-.00737917334364449*ii, -.023169798965672-.0115370319573312*ii, -.019796879433338-.0281761663938097*ii, -.0158990300749576-.0474048609604069*ii}
tw4 = {-.00609686563829517-.038745087669804*ii, .0104834650361434-.0240971080677847*ii, .0230950485407732-.0129553388502525*ii, .0236653481421159-.0124515046932896*ii, .0377594464968551+2.77555756156289e-17*ii, .048019227638035+.00906405715461864*ii, .0490033288920667+.00993346653975239*ii}
tw5 = {.0491685187212788+7.66747776381749e-16*ii, .0495940983990649-.0255916572433406*ii, .0498096050063981-.0385508550579585*ii, .0499649775475424-.0478939722785279*ii, .05-.05*ii}
twList2={tw1,tw2,tw3,tw4,tw5}
matLabTriangle2(2,last SetEncirclingTriangles)
apply(#twList2,i->matLabTL2(i+1,2,twList2_i))



{*
(Twist locus on a segment, (.05-.05*ii, .147249586032148-.047404860960407*ii), {.139278220184272-.0476175796315843*ii, .147249586032148-.047404860960407*ii})
tw1 = {.139278220184272-.0476175796315843*ii, .147249586032148-.047404860960407*ii}
(Twist locus on a segment, (.147249586032148-.047404860960407*ii, .212151944999172+.0099334665397523*ii), {.153497010396114-.0418855413252677*ii, .157739843322658-.0381371884505414*ii, .175200289648387-.0227116660846958*ii, .182992922881545-.0158272236157688*ii, .200561373056886-.000306284688375305*ii, .200660054056451-.000219104444680547*ii, .20066448058917-.000215193801265938*ii, .200717611456599-.000168255059354951*ii, .200908062603961+9.57567358739198e-16*ii, .212151944999172+.0099334665397523*ii})
tw2 = {.153497010396114-.0418855413252677*ii, .157739843322658-.0381371884505414*ii, .175200289648387-.0227116660846958*ii, .182992922881545-.0158272236157688*ii, .200561373056886-.000306284688375305*ii, .200660054056451-.000219104444680547*ii, .20066448058917-.000215193801265938*ii, .200717611456599-.000168255059354951*ii, .200908062603961+9.57567358739198e-16*ii, .212151944999172+.0099334665397523*ii}
(Twist locus on a segment, (.212151944999172+.0099334665397523*ii, .13004431729001+.0374713944206524*ii), {.203579515259291+.0128085581545913*ii, .199192077034142+.0142800530789608*ii, .198863231264089+.0143903440635172*ii, .192433236038032+.0165468884644459*ii, .13004431729001+.0374713944206524*ii})
tw3 = {.203579515259291+.0128085581545913*ii, .199192077034142+.0142800530789608*ii, .198863231264089+.0143903440635172*ii, .192433236038032+.0165468884644459*ii, .13004431729001+.0374713944206524*ii}
(Twist locus on a segment, (.13004431729001+.0374713944206524*ii, .147249586032148-.047404860960407*ii), {.133964753976512+.018131273066208*ii, .13764014498358-1.80411241501588e-16*ii, .140548235762022-.0143460622023175*ii, .144292991207221-.0328195207847156*ii, .146420755919154-.0433161136474567*ii, .147249586032148-.047404860960407*ii})
tw4 = {.133964753976512+.018131273066208*ii, .13764014498358-1.80411241501588e-16*ii, .140548235762022-.0143460622023175*ii, .144292991207221-.0328195207847156*ii, .146420755919154-.0433161136474567*ii, .147249586032148-.047404860960407*ii}
(Twist locus on a segment, (.147249586032148-.047404860960407*ii, .05-.05*ii), {.139372777271377-.0476150563427887*ii, .0620967876846337-.0496771930117645*ii, .05-.05*ii})
tw5 = {.139372777271377-.0476150563427887*ii, .0620967876846337-.0496771930117645*ii, .05-.05*ii}
Twist loci found! 
(END SEGMENT 0, .139278220184272-.0476175796315843*ii)
(END SEGMENT 1, .147249586032148-.047404860960407*ii)
The twist is : [5,6]
1 twists on fiber over: 
.1535-.041886*ii
{-1.0354+.07451*ii, .30881-1.5828*ii, .56404+2.1893*ii, .92294+1.3838*ii, 1.055-3.1389*ii, 1.055-2.0447*ii, 1.1295+3.1188*ii}
{-1.0354+.07451*ii, .30881-1.5828*ii, .56404+2.1893*ii, .92294+1.3838*ii, 1.055-2.0447*ii, 1.055-3.1389*ii, 1.1295+3.1188*ii}
(END SEGMENT 2, .153497010396114-.0418855413252677*ii)
(END SEGMENT 3, .157739843322658-.0381371884505414*ii)
(END SEGMENT 4, .175200289648387-.0227116660846958*ii)
(END SEGMENT 5, .182992922881545-.0158272236157688*ii)
The twist is : [6,5]
1 twists on fiber over: 
.20056-.000306*ii
{-1.1038+.00045*ii, .34816-1.9424*ii, .3494+1.9451*ii, 1.1011+1.7123*ii, 1.1016-1.7153*ii, 1.1016-3.1644*ii, 1.102+3.1642*ii}
{-1.1038+.00045*ii, .34816-1.9424*ii, .3494+1.9451*ii, 1.1011+1.7123*ii, 1.1016-3.1644*ii, 1.1016-1.7153*ii, 1.102+3.1642*ii}
(END SEGMENT 6, .200561373056886-.000306284688375305*ii)
The twist is : [6,7]
1 twists on fiber over: 
.20066-.000219*ii
{-1.104+.00032*ii, .3479-1.943*ii, .34879+1.9449*ii, 1.1016+1.7126*ii, 1.1017-3.1644*ii, 1.102-1.7148*ii, 1.102+3.1643*ii}
{-1.104+.00032*ii, .3479-1.943*ii, .34879+1.9449*ii, 1.1016+1.7126*ii, 1.1017-3.1644*ii, 1.102+3.1643*ii, 1.102-1.7148*ii}
(END SEGMENT 7, .200660054056451-.000219104444680547*ii)
The twist is : [5,4]
1 twists on fiber over: 
.20066-.000215*ii
{-1.104+.00031*ii, .34789-1.943*ii, .34876+1.9449*ii, 1.1017+1.7126*ii, 1.1017-3.1644*ii, 1.102+3.1643*ii, 1.102-1.7148*ii}
{-1.104+.00031*ii, .34789-1.943*ii, .34876+1.9449*ii, 1.1017-3.1644*ii, 1.1017+1.7126*ii, 1.102+3.1643*ii, 1.102-1.7148*ii}
(END SEGMENT 8, .20066448058917-.000215193801265938*ii)
The twist is : [5,6]
1 twists on fiber over: 
.20072-.000168*ii
{-1.1041+.00025*ii, .34775-1.9433*ii, .34843+1.9448*ii, 1.1017-3.1644*ii, 1.102+1.7128*ii, 1.102+3.1644*ii, 1.1022-1.7145*ii}
{-1.1041+.00025*ii, .34775-1.9433*ii, .34843+1.9448*ii, 1.1017-3.1644*ii, 1.102+3.1644*ii, 1.102+1.7128*ii, 1.1022-1.7145*ii}
(END SEGMENT 9, .200717611456599-.000168255059354951*ii)
The twist is : [2,3]
The twist is : [4,5]
The twist is : [7,6]
3 twists on fiber over: 
.20091
{-1.1044, .34725-1.9444*ii, .34725+1.9444*ii, 1.1019-3.1646*ii, 1.1019+3.1646*ii, 1.103+1.7135*ii, 1.103-1.7135*ii}
{-1.1044, .34725+1.9444*ii, .34725-1.9444*ii, 1.1019+3.1646*ii, 1.1019-3.1646*ii, 1.103-1.7135*ii, 1.103+1.7135*ii}
(END SEGMENT 10, .200908062603961+9.57567358739198e-16*ii)
(END SEGMENT 11, .212151944999172+.0099334665397523*ii)
The twist is : [5,6]
1 twists on fiber over: 
.20358+.012809*ii
{-1.1087-.01851*ii, .30576+1.8983*ii, .35529-2.0048*ii, 1.0929+3.1693*ii, 1.112-3.1642*ii, 1.112-1.65*ii, 1.1307+1.7699*ii}
{-1.1087-.01851*ii, .30576+1.8983*ii, .35529-2.0048*ii, 1.0929+3.1693*ii, 1.112-1.65*ii, 1.112-3.1642*ii, 1.1307+1.7699*ii}
(END SEGMENT 12, .203579515259291+.0128085581545913*ii)
The twist is : [5,4]
1 twists on fiber over: 
.19919+.01428*ii
{-1.1024-.02099*ii, .31838+1.8803*ii, .37635-2.0069*ii, 1.0908+3.1664*ii, 1.0908-1.6419*ii, 1.1124-3.1606*ii, 1.1137+1.7836*ii}
{-1.1024-.02099*ii, .31838+1.8803*ii, .37635-2.0069*ii, 1.0908-1.6419*ii, 1.0908+3.1664*ii, 1.1124-3.1606*ii, 1.1137+1.7836*ii}
(END SEGMENT 13, .199192077034142+.0142800530789608*ii)
The twist is : [6,7]
1 twists on fiber over: 
.19886+.01439*ii
{-1.102-.02118*ii, .31933+1.8789*ii, .37796-2.007*ii, 1.0891-1.6412*ii, 1.0906+3.1662*ii, 1.1125-3.1603*ii, 1.1125+1.7847*ii}
{-1.102-.02118*ii, .31933+1.8789*ii, .37796-2.007*ii, 1.0891-1.6412*ii, 1.0906+3.1662*ii, 1.1125+1.7847*ii, 1.1125-3.1603*ii}
(END SEGMENT 14, .198863231264089+.0143903440635172*ii)
The twist is : [6,5]
1 twists on fiber over: 
.19243+.016547*ii
{-1.0926-.02498*ii, .33784+1.8493*ii, .41074-2.0124*ii, 1.0561-1.6267*ii, 1.0874+3.1619*ii, 1.0874+1.8079*ii, 1.1131-3.155*ii}
{-1.0926-.02498*ii, .33784+1.8493*ii, .41074-2.0124*ii, 1.0561-1.6267*ii, 1.0874+1.8079*ii, 1.0874+3.1619*ii, 1.1131-3.155*ii}
(END SEGMENT 15, .192433236038032+.0165468884644459*ii)
(END SEGMENT 16, .13004431729001+.0374713944206524*ii)
The twist is : [3,4]
1 twists on fiber over: 
.13396+.018131*ii
{-.99215-.03636*ii, .46196+1.4564*ii, .71964-2.1801*ii, .71964-1.3715*ii, .91761+2.122*ii, 1.0689+3.1128*ii, 1.1044-3.1032*ii}
{-.99215-.03636*ii, .46196+1.4564*ii, .71964-1.3715*ii, .71964-2.1801*ii, .91761+2.122*ii, 1.0689+3.1128*ii, 1.1044-3.1032*ii}
(END SEGMENT 17, .133964753976512+.018131273066208*ii)
The twist is : [3,2]
The twist is : [4,5]
The twist is : [7,6]
3 twists on fiber over: 
.13764
{-.99755, .60321+1.4535*ii, .60321-1.4535*ii, .80816-2.116*ii, .80816+2.116*ii, 1.0874+3.1108*ii, 1.0874-3.1108*ii}
{-.99755, .60321-1.4535*ii, .60321+1.4535*ii, .80816+2.116*ii, .80816-2.116*ii, 1.0874-3.1108*ii, 1.0874+3.1108*ii}
(END SEGMENT 18, .13764014498358-1.80411241501588e-16*ii)
The twist is : [3,4]
1 twists on fiber over: 
.14055-.014346*ii
{-1.0043+.02775*ii, .4961-1.5008*ii, .71836+1.4283*ii, .71836+2.1365*ii, .89496-2.0844*ii, 1.0746-3.1177*ii, 1.1019+3.1103*ii}
{-1.0043+.02775*ii, .4961-1.5008*ii, .71836+2.1365*ii, .71836+1.4283*ii, .89496-2.0844*ii, 1.0746-3.1177*ii, 1.1019+3.1103*ii}
(END SEGMENT 19, .140548235762022-.0143460622023175*ii)
(END SEGMENT 20, .144292991207221-.0328195207847156*ii)
The twist is : [6,5]
1 twists on fiber over: 
.14642-.043316*ii
{-1.0236+.07975*ii, .3043-1.5424*ii, .58474+2.2143*ii, .90217+1.3471*ii, 1.051-2.0774*ii, 1.051-3.1336*ii, 1.1303+3.1122*ii}
{-1.0236+.07975*ii, .3043-1.5424*ii, .58474+2.2143*ii, .90217+1.3471*ii, 1.051-3.1336*ii, 1.051-2.0774*ii, 1.1303+3.1122*ii}
(END SEGMENT 21, .146420755919154-.0433161136474567*ii)
(END SEGMENT 22, .147249586032148-.047404860960407*ii)
(END SEGMENT 23, .139372777271377-.0476150563427887*ii)
(END SEGMENT 24, .0620967876846337-.0496771930117645*ii)
(END SEGMENT 25, .05-.05*ii)
END BRANCH POINT 


(Twist locus on a segment, (.05-.05*ii, .0490033288920667+.00993346653975239*ii), {.0499649775475425-.0478939722785278*ii, .0495940983990646-.0255916572433402*ii, .0491878892647111-.00116482138012693*ii, .0491685187212788+1.52655665885959e-16*ii, .0490136282555203+.00931412827396011*ii, .0490033288920667+.00993346653975239*ii})
tw1 = {.0499649775475425-.0478939722785278*ii, .0495940983990646-.0255916572433402*ii, .0491878892647111-.00116482138012693*ii, .0491685187212788+1.52655665885959e-16*ii, .0490136282555203+.00931412827396011*ii, .0490033288920667+.00993346653975239*ii}
(Twist locus on a segment, (.0490033288920667+.00993346653975239*ii, -.0331042988170953+.0374713944206525*ii), {.0297825501079181+.0163798885682567*ii, .0214277296070334+.019181996616556*ii, .0169808654119673+.0206734223086301*ii, .00515923208561742+.0246382584033228*ii, .00283970560009898+.02541620017125*ii, .00201954931854355+.0256912708940327*ii, -.00631591506102611+.0284868871269537*ii, -.0119704556044553+.0303833531939715*ii, -.0160530497270709+.0317526070092857*ii, -.0256515708655207+.0349718375515459*ii, -.0263694747059806+.0352126140145329*ii, -.0331042988170953+.0374713944206525*ii})
tw2 = {.0297825501079181+.0163798885682567*ii, .0214277296070334+.019181996616556*ii, .0169808654119673+.0206734223086301*ii, .00515923208561742+.0246382584033228*ii, .00283970560009898+.02541620017125*ii, .00201954931854355+.0256912708940327*ii, -.00631591506102611+.0284868871269537*ii, -.0119704556044553+.0303833531939715*ii, -.0160530497270709+.0317526070092857*ii, -.0256515708655207+.0349718375515459*ii, -.0263694747059806+.0352126140145329*ii, -.0331042988170953+.0374713944206525*ii}
(Twist locus on a segment, (-.0331042988170953+.0374713944206525*ii, -.0158990300749576-.0474048609604069*ii), {-.0306720609536044+.0254727883457755*ii, -.0306588803499729+.0254077663867164*ii, -.0286613257154425+.0155535200021312*ii, -.0285600422432051+.015053872947246*ii, -.0255084711235156-1.76941794549634e-15*ii, -.0240126386326793-.00737917334364449*ii, -.023169798965672-.0115370319573312*ii, -.019796879433338-.0281761663938097*ii, -.0158990300749576-.0474048609604069*ii})
tw3 = {-.0306720609536044+.0254727883457755*ii, -.0306588803499729+.0254077663867164*ii, -.0286613257154425+.0155535200021312*ii, -.0285600422432051+.015053872947246*ii, -.0255084711235156-1.76941794549634e-15*ii, -.0240126386326793-.00737917334364449*ii, -.023169798965672-.0115370319573312*ii, -.019796879433338-.0281761663938097*ii, -.0158990300749576-.0474048609604069*ii}
(Twist locus on a segment, (-.0158990300749576-.0474048609604069*ii, .0490033288920667+.00993346653975239*ii), {-.00609686563829517-.038745087669804*ii, .0104834650361434-.0240971080677847*ii, .0230950485407732-.0129553388502525*ii, .0236653481421159-.0124515046932896*ii, .0377594464968551+2.77555756156289e-17*ii, .048019227638035+.00906405715461864*ii, .0490033288920667+.00993346653975239*ii})
tw4 = {-.00609686563829517-.038745087669804*ii, .0104834650361434-.0240971080677847*ii, .0230950485407732-.0129553388502525*ii, .0236653481421159-.0124515046932896*ii, .0377594464968551+2.77555756156289e-17*ii, .048019227638035+.00906405715461864*ii, .0490033288920667+.00993346653975239*ii}
(Twist locus on a segment, (.0490033288920667+.00993346653975239*ii, .05-.05*ii), {.0491685187212788+7.66747776381749e-16*ii, .0495940983990649-.0255916572433406*ii, .0498096050063981-.0385508550579585*ii, .0499649775475424-.0478939722785279*ii, .05-.05*ii})
tw5 = {.0491685187212788+7.66747776381749e-16*ii, .0495940983990649-.0255916572433406*ii, .0498096050063981-.0385508550579585*ii, .0499649775475424-.0478939722785279*ii, .05-.05*ii}
Twist loci found! 
The twist is : [4,3]
1 twists on fiber over: 
.049965-.047894*ii
{-.80868+.17743*ii, .11182-1.0464*ii, .74187+2.5338*ii, .74187+.81338*ii, .9841-3.0396*ii, 1.0882-2.4357*ii, 1.1408+2.9971*ii}
{-.80868+.17743*ii, .11182-1.0464*ii, .74187+.81338*ii, .74187+2.5338*ii, .9841-3.0396*ii, 1.0882-2.4357*ii, 1.1408+2.9971*ii}
(END SEGMENT 0, .0499649775475425-.0478939722785278*ii)
The twist is : [5,6]
1 twists on fiber over: 
.049594-.025592*ii
{-.76961+.10533*ii, .22915-.96163*ii, .59037+.82605*ii, .81963+2.5182*ii, 1.0133-3.0164*ii, 1.0133-2.4633*ii, 1.1038+2.9919*ii}
{-.76961+.10533*ii, .22915-.96163*ii, .59037+.82605*ii, .81963+2.5182*ii, 1.0133-2.4633*ii, 1.0133-3.0164*ii, 1.1038+2.9919*ii}
(END SEGMENT 1, .0495940983990646-.0255916572433402*ii)
(END SEGMENT 2, .0491878892647111-.00116482138012693*ii)
The twist is : [2,3]
The twist is : [5,4]
The twist is : [6,7]
3 twists on fiber over: 
.049169
{-.74794, .39979-.87249*ii, .39979+.87249*ii, .91769+2.4955*ii, .91769-2.4955*ii, 1.0565-2.9965*ii, 1.0565+2.9965*ii}
{-.74794, .39979+.87249*ii, .39979-.87249*ii, .91769-2.4955*ii, .91769+2.4955*ii, 1.0565+2.9965*ii, 1.0565-2.9965*ii}
(END SEGMENT 3, .0491685187212788+1.52655665885959e-16*ii)
(END SEGMENT 4, .0490136282555203+.00931412827396011*ii)
(END SEGMENT 5, .0490033288920667+.00993346653975239*ii)
The twist is : [5,6]
1 twists on fiber over: 
.029783+.01638*ii
{-.66592-.09852*ii, .20109+.79159*ii, .49388-.67557*ii, .8694-2.5887*ii, 1.0077+2.5498*ii, 1.0077+2.972*ii, 1.0862-2.9506*ii}
{-.66592-.09852*ii, .20109+.79159*ii, .49388-.67557*ii, .8694-2.5887*ii, 1.0077+2.972*ii, 1.0077+2.5498*ii, 1.0862-2.9506*ii}
(END SEGMENT 6, .0297825501079181+.0163798885682567*ii)
(END SEGMENT 7, .0214277296070334+.019181996616556*ii)
(END SEGMENT 8, .0169808654119673+.0206734223086301*ii)
(END SEGMENT 9, .00515923208561742+.0246382584033228*ii)
(END SEGMENT 10, .00283970560009898+.02541620017125*ii)
(END SEGMENT 11, .00201954931854355+.0256912708940327*ii)
(END SEGMENT 12, -.00631591506102611+.0284868871269537*ii)
(END SEGMENT 13, -.0119704556044553+.0303833531939715*ii)
The twist is : [4,3]
1 twists on fiber over: 
-.016053+.031753*ii
{-.57804-.39389*ii, -.21917+.73291*ii, .781-.31132*ii, .781-2.7593*ii, .88305+2.9438*ii, 1.1727+2.6505*ii, 1.1795-2.8627*ii}
{-.57804-.39389*ii, -.21917+.73291*ii, .781-2.7593*ii, .781-.31132*ii, .88305+2.9438*ii, 1.1727+2.6505*ii, 1.1795-2.8627*ii}
(END SEGMENT 14, -.0160530497270709+.0317526070092857*ii)
The twist is : [7,6]
1 twists on fiber over: 
-.025652+.034972*ii
{-.5971-.44594*ii, -.28091+.75795*ii, .75882-2.7805*ii, .85334-.2826*ii, .86005+2.9506*ii, 1.2029+2.6574*ii, 1.2029-2.8568*ii}
{-.5971-.44594*ii, -.28091+.75795*ii, .75882-2.7805*ii, .85334-.2826*ii, .86005+2.9506*ii, 1.2029-2.8568*ii, 1.2029+2.6574*ii}
(END SEGMENT 15, -.0256515708655207+.0349718375515459*ii)
The twist is : [4,5]
1 twists on fiber over: 
-.026369+.035213*ii
{-.59867-.44946*ii, -.28507+.75986*ii, .75723-2.782*ii, .85845-.28087*ii, .85845+2.9511*ii, 1.2046-2.8565*ii, 1.205+2.6578*ii}
{-.59867-.44946*ii, -.28507+.75986*ii, .75723-2.782*ii, .85845+2.9511*ii, .85845-.28087*ii, 1.2046-2.8565*ii, 1.205+2.6578*ii}
(END SEGMENT 16, -.0263694747059806+.0352126140145329*ii)
(END SEGMENT 17, -.0331042988170953+.0374713944206525*ii)
The twist is : [6,7]
1 twists on fiber over: 
-.030672+.025473*ii
{-.5521-.48126*ii, -.32269+.70635*ii, .77383-2.8073*ii, .8465+2.9294*ii, .84682-.20399*ii, 1.2038-2.8316*ii, 1.2038+2.6884*ii}
{-.5521-.48126*ii, -.32269+.70635*ii, .77383-2.8073*ii, .8465+2.9294*ii, .84682-.20399*ii, 1.2038+2.6884*ii, 1.2038-2.8316*ii}
(END SEGMENT 18, -.0306720609536044+.0254727883457755*ii)
The twist is : [5,4]
1 twists on fiber over: 
-.030659+.025408*ii
{-.55173-.48128*ii, -.32273+.70594*ii, .77401-2.8074*ii, .84651+2.9293*ii, .84651-.2036*ii, 1.2037+2.6885*ii, 1.2037-2.8314*ii}
{-.55173-.48128*ii, -.32273+.70594*ii, .77401-2.8074*ii, .84651-.2036*ii, .84651+2.9293*ii, 1.2037+2.6885*ii, 1.2037-2.8314*ii}
(END SEGMENT 19, -.0306588803499729+.0254077663867164*ii)
The twist is : [3,4]
1 twists on fiber over: 
-.028661+.015554*ii
{-.492-.48925*ii, -.3338+.63996*ii, .80011-2.822*ii, .80011-.1377*ii, .84656+2.9033*ii, 1.1887+2.7147*ii, 1.1904-2.809*ii}
{-.492-.48925*ii, -.3338+.63996*ii, .80011-.1377*ii, .80011-2.822*ii, .84656+2.9033*ii, 1.1887+2.7147*ii, 1.1904-2.809*ii}
(END SEGMENT 20, -.0286613257154425+.0155535200021312*ii)
(END SEGMENT 21, -.0285600422432051+.015053872947246*ii)
The twist is : [1,2]
The twist is : [4,5]
The twist is : [7,6]
3 twists on fiber over: 
-.025508
{-.38398-.53034*ii, -.38398+.53034*ii, .74528, .83863-2.8557*ii, .83863+2.8557*ii, 1.1727+2.7637*ii, 1.1727-2.7637*ii}
{-.38398+.53034*ii, -.38398-.53034*ii, .74528, .83863+2.8557*ii, .83863-2.8557*ii, 1.1727-2.7637*ii, 1.1727+2.7637*ii}
(END SEGMENT 22, -.0255084711235156-1.76941794549634e-15*ii)
(END SEGMENT 23, -.0240126386326793-.00737917334364449*ii)
(END SEGMENT 24, -.023169798965672-.0115370319573312*ii)
The twist is : [3,4]
1 twists on fiber over: 
-.019797-.028176*ii
{-.55692+.41646*ii, -.24758-.70917*ii, .78543+.26845*ii, .78543+2.7754*ii, .87406-2.9342*ii, 1.1772-2.6667*ii, 1.1824+2.8498*ii}
{-.55692+.41646*ii, -.24758-.70917*ii, .78543+2.7754*ii, .78543+.26845*ii, .87406-2.9342*ii, 1.1772-2.6667*ii, 1.1824+2.8498*ii}
(END SEGMENT 25, -.019796879433338-.0281761663938097*ii)
(END SEGMENT 26, -.0158990300749576-.0474048609604069*ii)
The twist is : [4,3]
1 twists on fiber over: 
-.0060969-.038745*ii
{-.62582+.34468*ii, -.15408-.78909*ii, .77157+2.7208*ii, .77157+.40939*ii, .90352-2.9647*ii, 1.1608-2.6123*ii, 1.1725+2.8913*ii}
{-.62582+.34468*ii, -.15408-.78909*ii, .77157+.40939*ii, .77157+2.7208*ii, .90352-2.9647*ii, 1.1608-2.6123*ii, 1.1725+2.8913*ii}
(END SEGMENT 27, -.00609686563829517-.038745087669804*ii)
(END SEGMENT 28, .0104834650361434-.0240971080677847*ii)
The twist is : [5,6]
1 twists on fiber over: 
.023095-.012955*ii
{-.61847+.09404*ii, .18816-.7203*ii, .45262+.61285*ii, .888+2.6157*ii, 1.0059-2.9527*ii, 1.0059-2.5827*ii, 1.0779+2.9331*ii}
{-.61847+.09404*ii, .18816-.7203*ii, .45262+.61285*ii, .888+2.6157*ii, 1.0059-2.5827*ii, 1.0059-2.9527*ii, 1.0779+2.9331*ii}
(END SEGMENT 29, .0230950485407732-.0129553388502525*ii)
(END SEGMENT 30, .0236653481421159-.0124515046932896*ii)
The twist is : [2,3]
The twist is : [5,4]
The twist is : [6,7]
3 twists on fiber over: 
.037759
{-.69301, .36575-.78389*ii, .36575+.78389*ii, .93046+2.5407*ii, .93046-2.5407*ii, 1.0503-2.9737*ii, 1.0503+2.9737*ii}
{-.69301, .36575+.78389*ii, .36575-.78389*ii, .93046-2.5407*ii, .93046+2.5407*ii, 1.0503+2.9737*ii, 1.0503-2.9737*ii}
(END SEGMENT 31, .0377594464968551+2.77555756156289e-17*ii)
(END SEGMENT 32, .048019227638035+.00906405715461864*ii)
(END SEGMENT 33, .0490033288920667+.00993346653975239*ii)
The twist is : [3,2]
The twist is : [4,5]
The twist is : [7,6]
3 twists on fiber over: 
.049169
{-.74794, .39979+.87249*ii, .39979-.87249*ii, .91769-2.4955*ii, .91769+2.4955*ii, 1.0565+2.9965*ii, 1.0565-2.9965*ii}
{-.74794, .39979-.87249*ii, .39979+.87249*ii, .91769+2.4955*ii, .91769-2.4955*ii, 1.0565-2.9965*ii, 1.0565+2.9965*ii}
(END SEGMENT 34, .0491685187212788+7.66747776381749e-16*ii)
The twist is : [6,5]
1 twists on fiber over: 
.049594-.025592*ii
{-.76961+.10533*ii, .22915-.96163*ii, .59037+.82605*ii, .81963+2.5182*ii, 1.0133-2.4633*ii, 1.0133-3.0164*ii, 1.1038+2.9919*ii}
{-.76961+.10533*ii, .22915-.96163*ii, .59037+.82605*ii, .81963+2.5182*ii, 1.0133-3.0164*ii, 1.0133-2.4633*ii, 1.1038+2.9919*ii}
(END SEGMENT 35, .0495940983990649-.0255916572433406*ii)
(END SEGMENT 36, .0498096050063981-.0385508550579585*ii)
The twist is : [3,4]
1 twists on fiber over: 
.049965-.047894*ii
{-.80868+.17743*ii, .11182-1.0464*ii, .74187+.81338*ii, .74187+2.5338*ii, .9841-3.0396*ii, 1.0882-2.4357*ii, 1.1408+2.9971*ii}
{-.80868+.17743*ii, .11182-1.0464*ii, .74187+2.5338*ii, .74187+.81338*ii, .9841-3.0396*ii, 1.0882-2.4357*ii, 1.1408+2.9971*ii}
(END SEGMENT 37, .0499649775475424-.0478939722785279*ii)
(END SEGMENT 38, .05-.05*ii)
END BRANCH POINT 



*}





p1={.1535-.041886*ii,
.20056-.000306*ii,
.20066-.000219*ii,
.20066-.000215*ii,
.20072-.000168*ii,
.20091,
.20358+.012809*ii,
.19919+.01428*ii,
.19886+.01439*ii,
.19243+.016547*ii,
.13396+.018131*ii,
.13764,
.14055-.014346*ii,
.14642-.043316*ii}











p2={.049965-.047894*ii,
.049594-.025592*ii,
.049169,
.029783+.01638*ii,
-.016053+.031753*ii,
-.025652+.034972*ii,
-.026369+.035213*ii,
-.030672+.025473*ii,
-.030659+.025408*ii,
-.028661+.015554*ii,
-.025508,
-.019797-.028176*ii,
-.0060969-.038745*ii,
.023095-.012955*ii,
.037759,
.049169,
.049594-.025592*ii,
.049965-.047894*ii}


(p1|p2)/realPart//print
(p1|p2)/imaginaryPart//print
