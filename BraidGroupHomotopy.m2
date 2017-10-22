
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

sortBraid=method(TypicalValue=>Thing,Options=>{    	
    })
sortBraid(List) := o ->(twistedSols)->(
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
print segmentNumber;
oneSegment:={thePath_segmentNumber,thePath_(segmentNumber+1)};
twistLocusSegment:={};
gam':=oneSegment_0;
gam'':=oneSegment_1;
fs:=sub(f,{z=>x+ii*y,	t=>(1-s)*gam'+s*gam''	});
print 3;
g1:=value replace("ii","0",     toString  (fs));
g1=1/sub(max((flatten entries ((coefficients g1)_1))),CC)*g1;
print 1;
g2:=value replace("ii","0",     toString  (ii*fs));
g2=1/sub(max((flatten entries ((coefficients g2)_1))),CC)*g2;
print 2;
fiberG:={g1,g2,sub(g1,{y=>yA}),sub(g2,{y=>yA})};
makeB'InputFile(StoreBraidGroupFiles,B'Configs=>{"MPTYPE"=>2},    
    AffVariableGroup=>{x,y,yA,s},
    B'Polynomials=>fiberG    );
runBertini(StoreBraidGroupFiles);
realSols:=importSolutionsFile(StoreBraidGroupFiles,NameSolutionsFile=>"real_finite_solutions");
realSCoords:=radicalList((realSols/last),SetSameSCoordinateTolerance);
print (#realSols);
print radicalList(sort (realSols/last),SetSameSCoordinateTolerance);
branchS:={};
twistLocusSegment={};
for i in realSCoords do(
  if  (realPart i)<1 and 0<(realPart  i) then branchS=branchS|{ i});
if #branchS>0 then branchS=sort radicalList(branchS,SetSameSCoordinateTolerance);
twistLocusSegment=twistLocusSegment|radicalList((for i in branchS list 	sub(sub((1-s)*gam'+s*gam'',{s=>i}),CC))|{oneSegment_1},SetSameSCoordinateTolerance);
oneTwistLocusPathAndLoop=oneTwistLocusPathAndLoop|twistLocusSegment;
print (#oneTwistLocusPathAndLoop,segmentNumber));
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
    writeStartFile(StoreBraidGroupFiles,
	sortBraid importSolutionsFile(StoreBraidGroupFiles,	    NameSolutionsFile=>"raw_solutions_"|smallSegmentIndex,OrderPaths=>true)	);
    s1:= flatten importSolutionsFile(StoreBraidGroupFiles,
	    NameSolutionsFile=>"raw_solutions_"|smallSegmentIndex,OrderPaths=>true);
    s2:= flatten sortBraid importSolutionsFile(StoreBraidGroupFiles,
	    NameSolutionsFile=>"raw_solutions_"|smallSegmentIndex,OrderPaths=>true);
    if  (#radicalList((flatten s1)/realPart,1e-10)==#(flatten s1))
    then print "Untwisted fiber over: "
    else print (toString (#(flatten s1)-(#radicalList((flatten s1)/realPart,1e-10)))|" twists on fiber over: ");
--    print (#radicalList((flatten s1)/realPart,1e-10),#(flatten s1));
    if #radicalList((flatten s1)/realPart,1e-10)>0--<#(flatten s1) 
    then (
      printingPrecision=5;
      print (criticalTwistPointsOneBranchPoint_smallSegmentIndex);
      print (toString s1);
      print (toString s2);
      printingPrecision=300);
    moveB'File(StoreBraidGroupFiles,"final_parameters","start_parameters");        
    print ("END SEGMENT "|smallSegmentIndex)	);
    print ("END BRANCH POINT "))      


{*
--Example 1. 
restart
installPackage"BraidGroupHomotopy"
--
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
SetEncirclingTriangles==computeEncirclingTriangles(null)
#computeEncirclingTriangles(null)==#SetBraidGroupBranchPoints
oneTriangle=first SetEncirclingTriangles
--computeTwistLocusOfTriangle(varList,f,oneTriangle)
computeTwists(varList,f,oneTriangle)

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



------------------------------------------------------------------------------------------------------------------------------------------------------
--Determine all twist loci 
------------------------------------------------------------------------------------------------------------------------------------------------------
allTwistLoci={}
for branchPointIndex to #allBranchPointsT-1 do(
    print "branchPointIndex";
    print branchPointIndex;
oneBranchT=allBranchPointsT_branchPointIndex;
print oneBranchT;
oneTwistLocusPathAndLoop={};
thePath=theTriangles_branchPointIndex|startBasePoint;
for segmentNumber to #thePath-2 do(
print segmentNumber;
oneSegment={thePath_segmentNumber,thePath_(segmentNumber+1)};
twistLocusSegment={};
gam'=oneSegment_0;
gam''=oneSegment_1;
fs=sub(f,{z=>x+ii*y,	t=>(1-s)*gam'+s*gam''	});
print 3;
g1=value replace("ii","0",     toString  (fs));
g1=1/sub(max((flatten entries ((coefficients g1)_1))),CC)*g1;
print 1;
g2=value replace("ii","0",     toString  (ii*fs));
g2=1/sub(max((flatten entries ((coefficients g2)_1))),CC)*g2;
print 2;
fiberG={g1,g2,sub(g1,{y=>yA}),sub(g2,{y=>yA})};
makeB'InputFile(theDir2,B'Configs=>{{MPTYPE,2}},    AffVariableGroup=>{x,y,yA,s},       B'Polynomials=>fiberG    );
runBertini(theDir2);
readFile(theDir2);
realSols=importSolutionsFile(theDir2,NameSolutionsFile=>"real_finite_solutions");
realSCoords=radicalList((realSols/last),1e-10);
print (#realSols);
print radicalList(sort (realSols/last),1e-10);
branchS={};
twistLocusSegment={};
for i in realSCoords do(
  if  (realPart i)<1 and 0<(realPart  i) then branchS=branchS|{ i});
if #branchS>0 then branchS=sort radicalList(branchS,1e-10);
twistLocusSegment=twistLocusSegment|radicalList((for i in branchS list 	sub(sub((1-s)*gam'+s*gam'',{s=>i}),CC))|{oneSegment_1},1e-10);
oneTwistLocusPathAndLoop=oneTwistLocusPathAndLoop|twistLocusSegment;
print (#oneTwistLocusPathAndLoop,segmentNumber));
--
allTwistLoci=append(allTwistLoci,oneTwistLocusPathAndLoop))
---------------------------------------------------------------------------------

---------------------------------------------------------------------------------
--Now we determine the twist.
---------------------------------------------------------------------------------
writeStartFile(theDir, for i in startFiber list i,NameStartFile=>"start");
writeParameterFile(theDir, startBasePoint,NameParameterFile=>"start_parameters")
---------------
for branchPointIndex to #allTwistLoci-1 do(
criticalTwistPointsOneBranchPoint=allTwistLoci_branchPointIndex;--DOnt do radical List
for smallSegmentIndex to #criticalTwistPointsOneBranchPoint-1 do (
    writeParameterFile(theDir,
	{criticalTwistPointsOneBranchPoint_smallSegmentIndex},
	NameParameterFile=>"final_parameters");
    runBertini(theDir);
    moveB'File(theDir,"raw_solutions","raw_solutions_"|smallSegmentIndex);
    writeStartFile(theDir,
	sortBraid importSolutionsFile(theDir,	    NameSolutionsFile=>"raw_solutions_"|smallSegmentIndex,OrderPaths=>true)	);
    s1= flatten importSolutionsFile(theDir,
	    NameSolutionsFile=>"raw_solutions_"|smallSegmentIndex,OrderPaths=>true);
    s2= flatten sortBraid importSolutionsFile(theDir,
	    NameSolutionsFile=>"raw_solutions_"|smallSegmentIndex,OrderPaths=>true);
    if  (#radicalList((flatten s1)/realPart,1e-10)==#(flatten s1))
    then print "Untwisted fiber over: "
    else print (toString (#(flatten s1)-(#radicalList((flatten s1)/realPart,1e-10)))|" twists on fiber over: ");
--    print (#radicalList((flatten s1)/realPart,1e-10),#(flatten s1));
    if #radicalList((flatten s1)/realPart,1e-10)>0--<#(flatten s1) 
    then (
      printingPrecision=5;
      print (criticalTwistPointsOneBranchPoint_smallSegmentIndex);
      print (toString s1);
      print (toString s2);
      printingPrecision=300);
    moveB'File(theDir,"final_parameters","start_parameters");        
    print ("END SEGMENT "|smallSegmentIndex)	);
    print ("END BRANCH POINT "|branchPointIndex))      





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


