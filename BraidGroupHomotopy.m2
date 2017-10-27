
newPackage(
    "BraidGroupHomotopy",
    Version => "0.1", 
    Date => "October 15 2017",
    Authors => {
   {Name => "Jose Israel Rodriguez",
       Email => "JoIsRo@UChicago.edu",
       HomePage => "http://home.uchicago.edu/joisro/"}
   },
    Headline => "Using homotopy continuation to track braids. ",
    Configuration => { "BERTINIexecutable"=>"bertini" },
    DebuggingMode => true,
    AuxiliaryFiles => false,
    PackageImports => {"Bertini"},
    PackageExports => {}
)

  

exportMutable{
  "RotateTriangle",
  "RadiusBranchPoint",
  "StoreBraidGroupFiles",
  "SetSameXCoordinateTolerance",
  "SetSameTCoordinateTolerance",
  "SetSameSCoordinateTolerance",
  "SetBraidGroupBranchPoints",
  "SetDownstairsStartPoint",
  "SetUpstairsStartFiber",
  "SetEncirclingTriangles"
  }
RadiusBranchPoint=.05;

RotateTriangle=exp(.2*ii);
StoreBraidGroupFiles = temporaryFileName();
  makeDirectory StoreBraidGroupFiles 
SetSameXCoordinateTolerance=1e-8
SetSameTCoordinateTolerance=1e-8
SetSameSCoordinateTolerance=1e-8
SetBraidGroupBranchPoints={}
SetDownstairsStartPoint={}
SetUpstairsStartFiber={}



export {
--Methods
    "sortBraid",
    "computeBranchPoints",
    "computeBasePointAndFiber",
    "computeEncirclingTriangles",
    "computeTwistLocusOfTriangle",
    "computeBraid"
    }
 

--##########################################################################--
-- GLOBAL VARIABLES 
--##########################################################################--

BERTINIexe=(options BraidGroupHomotopy).Configuration#"BERTINIexecutable"
needsPackage "SimpleDoc"
printingPrecision=300
computeBranchPoints=method(TypicalValue=>Thing,Options=>{    	
    })
computeBranchPoints(Thing,Thing,Thing) := o ->(z,t,f)->(
    makeB'InputFile(StoreBraidGroupFiles,
    	AffVariableGroup=>{z,t},
	B'Configs=>{"MPTYPE"=>2,"SecurityLevel"=>1},
    	B'Polynomials=>{f,diff(z,f)}   );
    runBertini(StoreBraidGroupFiles);
    SetBraidGroupBranchPoints=radicalList(
	(importSolutionsFile(StoreBraidGroupFiles,NameSolutionsFile=>"finite_solutions"))/last//flatten,
	SetSameTCoordinateTolerance
	);
    return SetBraidGroupBranchPoints
    )
--We do not compute the set of nonpoperness. 

sortBraidOld=method(TypicalValue=>Thing,Options=>{    	
    })
sortBraidOld(List) := o ->(twistedSols)->(
    skipOne:=0;
    reorderedFiber:={};
    for onePoint to #twistedSols-1 do(	
      if onePoint<#twistedSols-1 then(
    	if abs((realPart first (twistedSols_onePoint))-(realPart first (twistedSols_(onePoint+1))))<SetSameXCoordinateTolerance 
    	then (
      	  reorderedFiber=reorderedFiber|{twistedSols_(onePoint+1),twistedSols_onePoint};
	  skipOne=1)
        else if skipOne===1 then skipOne=0 
	else reorderedFiber=reorderedFiber|{twistedSols_onePoint});
    if onePoint===#twistedSols-1 and skipOne===0 then reorderedFiber=(reorderedFiber|{twistedSols_onePoint})
    );
--    print reorderedFiber;
    if #reorderedFiber=!=#twistedSols then error"Bad Sort";
    return reorderedFiber)


sortBraid=method(TypicalValue=>Thing,Options=>{    	
    })
sortBraid(List) := o ->(twistedSols)->(
    twistedSols=twistedSols/first;
    reorderSols:={};
    solCounter:=0;
    while #twistedSols>0 do(
      if #twistedSols==1 
      then (	  
	  reorderSols=reorderSols|twistedSols;
	  twistedSols=drop(twistedSols,1))
      else 
        if #twistedSols>1 
	then(
          firstSol:=twistedSols_0;
	  secondSol:=twistedSols_1;
      	  if abs((realPart firstSol)-(realPart secondSol))>SetSameXCoordinateTolerance 
      	  then (
	    solCounter=solCounter+1;
	    twistedSols=drop(twistedSols,1);
	    reorderSols=reorderSols|{firstSol})
          else (
      	    if (imaginaryPart firstSol)<(imaginaryPart secondSol)
	    then   	print ("The twist is : "|"["|solCounter+1|","|solCounter+2|"]")
	    else   	print ("The twist is : "|"["|solCounter+2|","|solCounter+1|"]");
	    solCounter=solCounter+2;
	    twistedSols=drop(twistedSols,2);
	    reorderSols=reorderSols|{secondSol,firstSol})    		  
	  ));
      return apply(reorderSols,i->{i}))

---------------------------------------------------------------------------
computeBasePointAndFiber=method(TypicalValue=>Thing,Options=>{    	
    })
computeBasePointAndFiber(Thing,Thing,Thing) := o ->(z,t,f)->(
    makeB'InputFile(StoreBraidGroupFiles,
    B'Configs=>{
	"MPTYPE"=>2,
	"ParameterHomotopy"=>1},
    AffVariableGroup=>{z},
    ParameterGroup=>{t},
    B'Polynomials=>{f});
runBertini(StoreBraidGroupFiles,PreparePH2=>true);
SetUpstairsStartFiber= importSolutionsFile(StoreBraidGroupFiles)/last//sort;
writeStartFile(StoreBraidGroupFiles,apply(SetUpstairsStartFiber,i->{i}));
SetDownstairsStartPoint=importParameterFile(StoreBraidGroupFiles,NameParameterFile=>"start_parameters");
return(SetDownstairsStartPoint,SetUpstairsStartFiber)
  )
---------------------------------------------------------------------------
computeEncirclingTriangles=method(TypicalValue=>Thing,Options=>{    	
    })
computeEncirclingTriangles(Thing) := o ->(a)->(
SetEncirclingTriangles={};
for oneBranchT in SetBraidGroupBranchPoints do(
  topCorner:=oneBranchT+RadiusBranchPoint*exp(2*pi*ii/3)*RotateTriangle;
  leftCorner:=oneBranchT +RadiusBranchPoint*exp(4*pi*ii/3)*RotateTriangle;
  rightCorner:=oneBranchT+RadiusBranchPoint*RotateTriangle;
  if abs(leftCorner-SetDownstairsStartPoint_0)<abs(rightCorner-SetDownstairsStartPoint_0) and abs(leftCorner-SetDownstairsStartPoint_0)<abs(topCorner-SetDownstairsStartPoint_0) 
  then loopPoints:={leftCorner,rightCorner,topCorner,leftCorner};
  if abs(rightCorner-SetDownstairsStartPoint_0)<abs(leftCorner-SetDownstairsStartPoint_0) and abs(rightCorner-SetDownstairsStartPoint_0)<abs(topCorner-SetDownstairsStartPoint_0) 
  then loopPoints={rightCorner,topCorner,leftCorner,rightCorner};--
  if abs(topCorner-SetDownstairsStartPoint_0)<abs(leftCorner-SetDownstairsStartPoint_0) and abs(topCorner-SetDownstairsStartPoint_0)<abs(rightCorner-SetDownstairsStartPoint_0) 
  then loopPoints={topCorner,leftCorner,rightCorner,topCorner};
  SetEncirclingTriangles=append(SetEncirclingTriangles,loopPoints));
return SetEncirclingTriangles)


computeTwistLocusOfTriangle=method(TypicalValue=>Thing,Options=>{    	
    })
computeTwistLocusOfTriangle(Thing,Thing,Thing) := o ->(varList,f,aTriangle)->(
(z,t,x,y,yA,s):=varList;
oneTwistLocusPathAndLoop:={};
thePath:=aTriangle|SetDownstairsStartPoint;
for segmentNumber to #thePath-1-1 do(
--print segmentNumber;
oneSegment:={thePath_segmentNumber,thePath_(segmentNumber+1)};
twistLocusSegment:={};
gam':=oneSegment_0;
gam'':=oneSegment_1;
fs:=sub(f,{z=>x+ii*y,	t=>(1-s)*gam'+s*gam''	});
--print 3;
g1:=value replace("ii","0",     toString  (fs));
g1=1/sub(max((flatten entries ((coefficients g1)_1))),CC)*g1;
--print 1;
g2:=value replace("ii","0",     toString  (ii*fs));
g2=1/sub(max((flatten entries ((coefficients g2)_1))),CC)*g2;
--print 2;
fiberG:={g1,g2,sub(g1,{y=>yA}),sub(g2,{y=>yA})};
makeB'InputFile(StoreBraidGroupFiles,B'Configs=>{"MPTYPE"=>2},    
    AffVariableGroup=>{x,y,yA,s},
    B'Polynomials=>fiberG    );
runBertini(StoreBraidGroupFiles);
realSols:=importSolutionsFile(StoreBraidGroupFiles,NameSolutionsFile=>"real_finite_solutions");
realSCoords:=radicalList((realSols/last),SetSameSCoordinateTolerance);
--print (#realSols);
--print radicalList(sort (realSols/last),SetSameSCoordinateTolerance);
branchS:={};
twistLocusSegment={};
for i in realSCoords do(
  if  (realPart i)<1 and 0<(realPart  i) then branchS=branchS|{ i});
if #branchS>0 then branchS=sort radicalList(branchS,SetSameSCoordinateTolerance);
twistLocusSegment=twistLocusSegment|radicalList((for i in branchS list 	sub(sub((1-s)*gam'+s*gam'',{s=>i}),CC))|{oneSegment_1},SetSameSCoordinateTolerance);
oneTwistLocusPathAndLoop=oneTwistLocusPathAndLoop|twistLocusSegment;
--print (#oneTwistLocusPathAndLoop,segmentNumber)
);
print ("Twist loci found! ");
return oneTwistLocusPathAndLoop
)


--Now we determine the twist.
---------------------------------------------------------------------------------
computeBraid=method(TypicalValue=>Thing,Options=>{    	
    })
computeBraid(Thing,Thing,Thing) := o ->(varList,f,aTriangle)->(
    (z,t,x,y,yA,s):=varList;
    criticalTwistPointsOneBranchPoint:=computeTwistLocusOfTriangle(varList,f,aTriangle);
    writeStartFile(StoreBraidGroupFiles, for i in SetUpstairsStartFiber list {i} ,NameStartFile=>"start");
    writeParameterFile(StoreBraidGroupFiles, SetDownstairsStartPoint,NameParameterFile=>"start_parameters");    
    makeB'InputFile(StoreBraidGroupFiles,
    B'Configs=>{
	"MPTYPE"=>2,
	"MaxSecurityLevel"=>1,
	"ParameterHomotopy"=>2},
    AffVariableGroup=>{z},
    ParameterGroup=>{t},
    B'Polynomials=>{f});
for smallSegmentIndex to #criticalTwistPointsOneBranchPoint-1 do (
    writeParameterFile(StoreBraidGroupFiles,
	{criticalTwistPointsOneBranchPoint_smallSegmentIndex},
	NameParameterFile=>"final_parameters");
    runBertini(StoreBraidGroupFiles);
    moveB'File(StoreBraidGroupFiles,"raw_solutions","raw_solutions_"|smallSegmentIndex);
    s1Points:= importSolutionsFile(StoreBraidGroupFiles,
	    NameSolutionsFile=>"raw_solutions_"|smallSegmentIndex,OrderPaths=>true);
    s1:=flatten s1Points;
    s2Points:=sortBraid(s1Points);
    s2:=flatten s2Points;
    writeStartFile(StoreBraidGroupFiles,s2Points);
    if not (#radicalList((flatten s1)/realPart,SetSameXCoordinateTolerance)==#(flatten s1))
    then print (toString (#(flatten s1)-(#radicalList((flatten s1)/realPart,SetSameXCoordinateTolerance)))|" twists on fiber over: ");
--    print (#radicalList((flatten s1)/realPart,1e-10),#(flatten s1));
    if #radicalList((flatten s1)/realPart,1e-10)<#(flatten s1) ---<0 to print allfibers
    then (
      printingPrecision=5;
      print (criticalTwistPointsOneBranchPoint_smallSegmentIndex);
      print (toString s1);
      print (toString s2);
      printingPrecision=300);
    moveB'File(StoreBraidGroupFiles,"final_parameters","start_parameters");        
--    print ("END SEGMENT "|smallSegmentIndex)	
);
    print ("END BRANCH POINT "))      


{*
--Example 1. 
restart
installPackage"BraidGroupHomotopy"
--
printingPrecision=300
R=CC[z,t,x,y,yA,s]
varList=(z,t,x,y,yA,s)
f=z^3-t^2; 
computeBranchPoints(z,t,f)
SetBraidGroupBranchPoints
cbpaf=computeBasePointAndFiber(z,t,f)
SetDownstairsStartPoint==first cbpaf
SetUpstairsStartFiber==last cbpaf
#SetDownstairsStartPoint==1
#SetUpstairsStartFiber==3
SetEncirclingTriangles=computeEncirclingTriangles(null)
#computeEncirclingTriangles(null)==#SetBraidGroupBranchPoints
oneTriangle=first SetEncirclingTriangles
--computeTwistLocusOfTriangle(varList,f,oneTriangle)
computeBraid(varList,f,oneTriangle)


i32 : --computeTwistLocusOfTriangle(varList,f,oneTriangle)
      computeBraid(varList,f,oneTriangle)
Twist loci found! 
The twist is : [1,2]
1 twists on fiber over: 
-.025508
{-.043327-.075045*ii, -.043327+.075045*ii, .086654}
{-.043327-.075045*ii, -.043327+.075045*ii, .086654}
The twist is : [2,3]
1 twists on fiber over: 
-.033359*ii
{-.10363, .051814-.089744*ii, .051814+.089744*ii}
{-.10363, .051814-.089744*ii, .051814+.089744*ii}
The twist is : [1,2]
1 twists on fiber over: 
.037759
{-.056276-.097472*ii, -.056276+.097472*ii, .11255}
{-.056276-.097472*ii, -.056276+.097472*ii, .11255}
The twist is : [2,3]
1 twists on fiber over: 
.026369*ii
{-.088591, .044296-.076722*ii, .044296+.076722*ii}
{-.088591, .044296-.076722*ii, .044296+.076722*ii}
END BRANCH POINT 



--Example 2. 
restart
installPackage"BraidGroupHomotopy"
--
printingPrecision=300
R=CC[z,t,x,y,yA,s]
varList=(z,t,x,y,yA,s)
f=z^4-4*z^2+3+t; 
--allBranchPointsT={ -3,1};
computeBranchPoints(z,t,f)
SetBraidGroupBranchPoints
cbpaf=computeBasePointAndFiber(z,t,f)
SetDownstairsStartPoint==first cbpaf
SetUpstairsStartFiber==last cbpaf
#SetDownstairsStartPoint==1
#SetUpstairsStartFiber==first degree f
SetEncirclingTriangles=computeEncirclingTriangles(null)
#computeEncirclingTriangles(null)==#SetBraidGroupBranchPoints
oneTriangle=first SetEncirclingTriangles
--computeTwistLocusOfTriangle(varList,f,oneTriangle)
computeBraid(varList,f,oneTriangle)

twoTriangle=last SetEncirclingTriangles;
computeBraid(varList,f,twoTriangle)


i16 : --computeTwistLocusOfTriangle(varList,f,oneTriangle)
      computeBraid(varList,f,oneTriangle)
Twist loci found! 
The twist is : [1,2]
The twist is : [3,4]
2 twists on fiber over: 
1.0378
{-1.4159-.06862*ii, -1.4159+.06862*ii, 1.4159-.06862*ii, 1.4159+.06862*ii}
{-1.4159-.06862*ii, -1.4159+.06862*ii, 1.4159-.06862*ii, 1.4159+.06862*ii}
END BRANCH POINT 


ii18 : twoTriangle=last SetEncirclingTriangles;

ii19 : computeBraid(varList,f,twoTriangle)
Twist loci found! 
The twist is : [2,3]
1 twists on fiber over: 
-3.0255
{-2.0016, -.079793*ii, .079793*ii, 2.0016}
{-2.0016, -.079793*ii, .079793*ii, 2.0016}
END BRANCH POINT 

--Example 3. 
restart
installPackage"BraidGroupHomotopy"
--
printingPrecision=300
R=CC[z,t,x,y,yA,s]
varList=(z,t,x,y,yA,s)
f=z^6-4*z^3+3+t
computeBranchPoints(z,t,f)
SetBraidGroupBranchPoints
cbpaf=computeBasePointAndFiber(z,t,f)
SetDownstairsStartPoint==first cbpaf
SetUpstairsStartFiber==last cbpaf
#SetDownstairsStartPoint==1
#SetUpstairsStartFiber==first degree f
SetEncirclingTriangles=computeEncirclingTriangles(null)
#computeEncirclingTriangles(null)==#SetBraidGroupBranchPoints
oneTriangle=first SetEncirclingTriangles
--computeTwistLocusOfTriangle(varList,f,oneTriangle)
computeBraid(varList,f,oneTriangle)



---Example 4
--
restart
installPackage"BraidGroupHomotopy"
--
printingPrecision=300
R=CC[z,t,x,y,yA,s]
varList=(z,t,x,y,yA,s)
f=z^4-2*t^3*z^2-4*t^5*z+t^6-t^7
computeBranchPoints(z,t,f)
SetBraidGroupBranchPoints
cbpaf=computeBasePointAndFiber(z,t,f)
SetDownstairsStartPoint==first cbpaf
SetUpstairsStartFiber==last cbpaf
#SetDownstairsStartPoint==1
#SetUpstairsStartFiber==first degree f
SetEncirclingTriangles=computeEncirclingTriangles(null)
#computeEncirclingTriangles(null)==#SetBraidGroupBranchPoints
oneTriangle=first SetEncirclingTriangles
--computeTwistLocusOfTriangle(varList,f,oneTriangle)
computeBraid(varList,f,oneTriangle)




--Botong's emailed example:
f=z^4-2*z^3*t^2-4*z^5*t+z^6-z^7
f=sub(f,{z=>t,t=>z})
makeB'InputFile(theDir,
    AffVariableGroup=>{t,z},
    B'Polynomials=>{f,diff(z,f)}   )
runBertini(theDir)
allBranchPointsT=radicalList( (importSolutionsFile(theDir))/first,1e-10)

----Example 5 \label{ex:dessinA}
restart
installPackage"BraidGroupHomotopy"
--
printingPrecision=300
R=CC[z,t,x,y,yA,s]
varList=(z,t,x,y,yA,s)
alpha=(34+6*sqrt(21))/7
--alpha=(34-6*sqrt(21))/7
f=1/100*z^3*(z^2-2*z+alpha)^2-t
computeBranchPoints(z,t,f)--there should be two branch points. 
--{-1.6314861610711+8.10462807976364e-15*ii, -8.69326802036999e-14+3.23872870621213e-14*ii}
#SetBraidGroupBranchPoints===2

cbpaf=computeBasePointAndFiber(z,t,f)

SetDownstairsStartPoint==first cbpaf
SetUpstairsStartFiber==last cbpaf
#SetDownstairsStartPoint==1
#SetUpstairsStartFiber==first degree f
SetEncirclingTriangles=computeEncirclingTriangles(null)
#computeEncirclingTriangles(null)==#SetBraidGroupBranchPoints

theBraids=apply(SetEncirclingTriangles,oneTriangle->computeBraid(varList,f,oneTriangle))




----Example 5 \label{ex:dessinB}
restart
installPackage"BraidGroupHomotopy"
--
printingPrecision=300
R=CC[z,t,x,y,yA,s]
varList=(z,t,x,y,yA,s)
beta=(34-6*sqrt(21))/7
f=z^3*(z^2-2*z+beta)^2-t/20
computeBranchPoints(z,t,f)--there should be two branch points. 

#SetBraidGroupBranchPoints===2

cbpaf=computeBasePointAndFiber(z,t,f)

SetDownstairsStartPoint==first cbpaf
SetUpstairsStartFiber==last cbpaf
#SetDownstairsStartPoint==1
#SetUpstairsStartFiber==first degree f
SetEncirclingTriangles=computeEncirclingTriangles(null)
#computeEncirclingTriangles(null)==#SetBraidGroupBranchPoints

theBraids=apply(SetEncirclingTriangles,oneTriangle->computeBraid(varList,f,oneTriangle))


*}



///




--##########################################################################--
-- TESTS
--##########################################################################--

--TEST///
--load concatenate(Bertini#"source directory","./Bertini/TST/bertiniIsProjective.tst.m2")
--///



---newtst

--##########################################################################--
-- DOCUMENTATION
--##########################################################################--

beginDocumentation()

load "./BraidGroupHomotopySupplement/doc.m2";
end




{*
This is the non generic example
restart
R=CC[z,x,y,yA,t,s,a,da,b,db]
f=(z^3-t)*(z^3+1/8*t)
fPrime=diff(z,f)
fOnLine=sub(f,{t=>random RR+s*1+random RR*ii+s*ii})
fOnLine=sub(f,{t=>-.5+s-.25*ii})
xyfOnLine=sub(fOnLine,{z=>x+ii*y})
xyBarfOnLine=sub(fOnLine,{z=>x+ii*yA})
G1=value replace("ii","0", toString (xyfOnLine))
G2=value replace("ii","0", toString (ii*xyfOnLine))
G1Bar=value replace("ii","0", toString (xyBarfOnLine))
G2Bar=value replace("ii","0", toString (ii*xyBarfOnLine))
loadPackage"Bertini"
win=bertiniZeroDimSolve({G1,G2,G1Bar,G2Bar},AffVariableGroup=>{x,y,yA,s});
win=oo;
delete(null,apply(win/coordinates,i->if abs imaginaryPart(last i)<1e-5 and 0<(realPart(last i))and  1>(realPart(last i)) then i else null))
netList oo
*}










--restart
printingPrecision=300
R=CC[z,x,y,yA,t,s]
---
f=z^3-t^2; allBranchPointsT={ 0};
---
f=z^4-4*z^2+3+t;allBranchPointsT={ -3,1};
---
f=z^6-4*z^3+3+t
makeB'InputFile(theDir,
    AffVariableGroup=>{z},
    B'Polynomials=>{diff(z,f)}   )
runBertini(theDir)
allBranchPointsT=radicalList(flatten importSolutionsFile(theDir),1e-10)
--
f=z^4-2*t^3*z^2-4*t^5*z+t^6-t^7
makeB'InputFile(theDir,
    AffVariableGroup=>{t,z},
    B'Polynomials=>{f,diff(z,f)}   )
runBertini(theDir)
allBranchPointsT=radicalList( (importSolutionsFile(theDir))/first,1e-10)
---

--Botong's emailed example:
f=z^4-2*z^3*t^2-4*z^5*t+z^6-z^7
f=sub(f,{z=>t,t=>z})
makeB'InputFile(theDir,
    AffVariableGroup=>{t,z},
    B'Polynomials=>{f,diff(z,f)}   )
runBertini(theDir)
allBranchPointsT=radicalList( (importSolutionsFile(theDir))/first,1e-10)


L={-.32844+.055715*ii, .10536-.34454*ii, .2258+.2872*ii, .96477+2.7327*ii, 1.0007-2.7228*ii, 1.0007-2.8521*ii, 1.0311+2.8437*ii}
netList sortBraid(apply(L,i->{i})    )
