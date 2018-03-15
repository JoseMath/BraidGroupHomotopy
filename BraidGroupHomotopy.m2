
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
  "SetRotateTriangle",
  "RadiusBranchPoint",
  "StoreBraidGroupFiles",
  "SetSameXCoordinateTolerance",
  "SetSameTCoordinateTolerance",
  "SetSameSCoordinateTolerance",
  "SetBraidGroupBranchPoints",
  "SetDownstairsStartPoint",
  "SetUpstairsStartFiber",
  "SetEncirclingTriangles",
  "SetHyperplaneArrangement"
--Options
  }
RadiusBranchPoint=.05;

SetRotateTriangle=exp(.2*ii);
StoreBraidGroupFiles = temporaryFileName();
  makeDirectory StoreBraidGroupFiles 
SetSameXCoordinateTolerance=1e-8
SetSameTCoordinateTolerance=1e-8
SetSameSCoordinateTolerance=1e-8
SetBraidGroupBranchPoints={}
SetDownstairsStartPoint={}
SetUpstairsStartFiber={}
SetHyperplaneArrangement={}


export {
--Methods
    "sortBraid",
    "computeBranchPoints",
    "computeBasePointAndFiber",
    "computeEncirclingTriangles",
    "computeTwistLocusOfTriangle",
    "computeBraid",
    "IsHyperplaneArrangement",
    "LatexPrint"
    }
 

--##########################################################################--
-- GLOBAL VARIABLES 
--##########################################################################--

BERTINIexe=(options BraidGroupHomotopy).Configuration#"BERTINIexecutable"
needsPackage "SimpleDoc"
printingPrecision=300

---This function computes the branch points of a univariate polynomial in one variable. 
--z is the unknown and t is the parameter. 
--It changes the value of SetBraidGroupBranchPoints to a list of complex numbers.
--We do not compute the set of nonpoperness and assume genericness in this regard. 
computeBranchPoints=method(TypicalValue=>Thing,Options=>{ })
computeBranchPoints(Thing,Thing,Thing) := o ->(z,t,f)->(
    makeB'InputFile(StoreBraidGroupFiles,
    	AffVariableGroup=>{z,t},
	B'Configs=>{"MPTYPE"=>2,"SecurityLevel"=>1},
    	B'Polynomials=>{f,diff(z,f)}   );
    runBertini(StoreBraidGroupFiles);
    SetBraidGroupBranchPoints=radicalList(
	(importSolutionsFile(StoreBraidGroupFiles,NameSolutionsFile=>"finite_solutions")
	  )/last//flatten,SetSameTCoordinateTolerance	);
    return SetBraidGroupBranchPoints    )



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


sortBraid=method(TypicalValue=>Thing,Options=>{})
--Input a list of complex numbers or list of numbers that are nested in a list.
--Outputs a sorted list of complex numbers.
sortBraid(List) := o ->(twistedSols)->(
    twistedSols=flatten twistedSols;
    reorderSols:={};
    solCounter:=0;
    while #twistedSols>0 do(
      if #twistedSols==1 
      then (	  
	reorderSols=reorderSols|twistedSols;
	twistedSols=drop(twistedSols,1))
      else (
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
	    then (
		print ("The twist is : "|"["|solCounter+1|","|solCounter+2|"]");
		print("a_{"|solCounter+1|"} ") 		)
	    else (
		print ("The twist is : "|"["|solCounter+2|","|solCounter+1|"]");
		print("a_{"|solCounter+1|"}^{-1} ") 		);
	    solCounter=solCounter+2;
	    twistedSols=drop(twistedSols,2);
	    reorderSols=reorderSols|{secondSol,firstSol})    		  
	  )));
      return apply(reorderSols,i->{i}))

sortComplexNumbers=twistedSols->(
--Verbose    print twistedSols;
    twistedSols=sort twistedSols;
--Verbose    print twistedSols;
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
          firstSol :=twistedSols_0;
	  secondSol:=twistedSols_1;
      	  if abs((realPart firstSol)-(realPart secondSol))>SetSameXCoordinateTolerance 
      	  then (
	    solCounter =solCounter+1;
	    twistedSols=drop(twistedSols,1);
	    reorderSols=reorderSols|{firstSol})
          else (
      	    if (imaginaryPart firstSol)<(imaginaryPart secondSol)
	    then reorderSols=reorderSols|{firstSol,secondSol}
	    else reorderSols=reorderSols|{secondSol,firstSol};
	    solCounter=solCounter+2;
	    twistedSols=drop(twistedSols,2))));
      print reorderSols;
      return reorderSols)    
    
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
SetUpstairsStartFiber= sortComplexNumbers flatten (importSolutionsFile(StoreBraidGroupFiles));
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
  topCorner:=oneBranchT+RadiusBranchPoint*exp(2*pi*ii/3)*SetRotateTriangle;
  leftCorner:=oneBranchT +RadiusBranchPoint*exp(4*pi*ii/3)*SetRotateTriangle;
  rightCorner:=oneBranchT+RadiusBranchPoint*SetRotateTriangle;
  if abs(leftCorner-SetDownstairsStartPoint_0)<abs(rightCorner-SetDownstairsStartPoint_0) and abs(leftCorner-SetDownstairsStartPoint_0)<abs(topCorner-SetDownstairsStartPoint_0) 
  then loopPoints:={leftCorner,rightCorner,topCorner,leftCorner};
  if abs(rightCorner-SetDownstairsStartPoint_0)<abs(leftCorner-SetDownstairsStartPoint_0) and abs(rightCorner-SetDownstairsStartPoint_0)<abs(topCorner-SetDownstairsStartPoint_0) 
  then loopPoints={rightCorner,topCorner,leftCorner,rightCorner};--
  if abs(topCorner-SetDownstairsStartPoint_0)<abs(leftCorner-SetDownstairsStartPoint_0) and abs(topCorner-SetDownstairsStartPoint_0)<abs(rightCorner-SetDownstairsStartPoint_0) 
  then loopPoints={topCorner,leftCorner,rightCorner,topCorner};
  SetEncirclingTriangles=append(SetEncirclingTriangles,loopPoints));
return SetEncirclingTriangles)




computeTwistLocusOfTriangle=method(TypicalValue=>Thing,Options=>{ 
	IsHyperplaneArrangement=>false   	
    })
computeTwistLocusOfTriangle(Thing,Thing,Thing) := o ->(varList,f,aTriangle)->(
(z,t,x,y,yA,s):=varList;
oneTwistLocusPathAndLoop:={};
thePath:=SetDownstairsStartPoint|aTriangle|SetDownstairsStartPoint;
for segmentNumber to #thePath-1-1 do(
--print segmentNumber;
oneSegment:={thePath_segmentNumber,thePath_(segmentNumber+1)};
gam':=oneSegment_0;
gam'':=oneSegment_1;
if not o.IsHyperplaneArrangement 
then (
  fs:=sub(f,{z=>x+ii*y,	t=>(1-s)*gam'+s*gam''	});
  g1:=value replace("ii","0",     toString  (fs));
  g2:=value replace("ii","0",     toString  (ii*fs));
  fiberG:={g1,g2,sub(g1,{y=>yA}),sub(g2,{y=>yA})};
  makeB'InputFile(StoreBraidGroupFiles,B'Configs=>{"MPTYPE"=>2},    
    AffVariableGroup=>{x,y,yA,s},
    B'Polynomials=>fiberG    );
  runBertini(StoreBraidGroupFiles);
  realSols:=importSolutionsFile(StoreBraidGroupFiles,NameSolutionsFile=>"real_finite_solutions");
  realSCoords:=radicalList((realSols/last),SetSameSCoordinateTolerance);
--print (#realSols);
--print radicalList(sort (realSols/last),SetSameSCoordinateTolerance);
  branchS:=delete(null,apply(realSCoords,i->if  (realPart i)<1 and 0<(realPart  i) then i else null));
  if #branchS>0 then branchS=sort branchS;
  twistLocusSegment:=    (for i in branchS list 	sub(sub((1-s)*gam'+s*gam'',{s=>i}),CC))|{gam''}    ;
)
else if o.IsHyperplaneArrangement then(
  branchS={};
--  print 1;
  apply(subsets(SetHyperplaneArrangement,2),H2->(
  fs=sub(product H2,{z=>x+ii*y,	t=>(1-s)*gam'+s*gam''	});
--  print fs;
  g1=value replace("ii","0",     toString  (fs));
  g2=value replace("ii","0",     toString  (ii*fs));
  fiberG={g1,g2,sub(g1,{y=>yA}),sub(g2,{y=>yA})};
  makeB'InputFile(StoreBraidGroupFiles,B'Configs=>{"MPTYPE"=>2},    
    AffVariableGroup=>{x,y,yA,s},
    B'Polynomials=>fiberG    );
  runBertini(StoreBraidGroupFiles);
  realSols=importSolutionsFile(StoreBraidGroupFiles,NameSolutionsFile=>"real_finite_solutions");
  realSCoords=radicalList((realSols/last),SetSameSCoordinateTolerance);
--print (#realSols);
--print radicalList(sort (realSols/last),SetSameSCoordinateTolerance);
  branchS=branchS|delete(null,apply(realSCoords,i->if  (realPart i)<1 and 0<(realPart  i) then i else null));
  if #branchS>0 then branchS=sort radicalList( branchS,SetSameSCoordinateTolerance);
  ));
  twistLocusSegment=    (for i in branchS list 	sub(sub((1-s)*gam'+s*gam'',{s=>i}),CC))|{gam''}    ;
    );
print ("Twist locus on a segment",(gam',gam''),twistLocusSegment);
print("tw"|segmentNumber+1|" = "|toString twistLocusSegment);
oneTwistLocusPathAndLoop=oneTwistLocusPathAndLoop|twistLocusSegment;
--print (#oneTwistLocusPathAndLoop,segmentNumber)
);
print ("Twist loci found! ");
return oneTwistLocusPathAndLoop
)



--Now we determine the twist.
---------------------------------------------------------------------------------
computeBraid=method(TypicalValue=>Thing,Options=>{    	
    IsHyperplaneArrangement=>false   	
    })
computeBraid(Thing,Thing,Thing) := o ->(varList,f,aTriangle)->(
    (z,t,x,y,yA,s):=varList;
    criticalTwistPointsOneBranchPoint:=computeTwistLocusOfTriangle(varList,f,aTriangle,IsHyperplaneArrangement=>o.IsHyperplaneArrangement);
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
  s1Points:= importSolutionsFile(StoreBraidGroupFiles,
	    NameSolutionsFile=>"raw_solutions",OrderPaths=>true);
  s1:=flatten s1Points;
  moveB'File(StoreBraidGroupFiles,"raw_solutions","raw_solutions_"|smallSegmentIndex);
  s2Points:=sortBraid(s1Points);
  s2:=flatten s2Points;
  writeStartFile(StoreBraidGroupFiles,s2Points);
  if not (#radicalList((flatten s1)/realPart,SetSameXCoordinateTolerance)==#(flatten s1))
  then print (toString (#(flatten s1)-(#radicalList((flatten s1)/realPart,SetSameXCoordinateTolerance)))|" twists on fiber over: ");
--    print (#radicalList((flatten s1)/realPart,1e-10),#(flatten s1));
  if #radicalList((flatten s1)/realPart,1e-10)<#(flatten s1) ---<0 to print allfibers
  then (
      printingPrecision=5;
--VERBOSE
--      print (criticalTwistPointsOneBranchPoint_smallSegmentIndex);
--      print (toString s1);
--      print (toString s2);
      printingPrecision=300);
  moveB'File(StoreBraidGroupFiles,"final_parameters","start_parameters");        
--VERBOSE  print ("END SEGMENT "|smallSegmentIndex,criticalTwistPointsOneBranchPoint_smallSegmentIndex)	
);
    print ("END BRANCH POINT "))      


---Now we have a new set of functions specifically designed for hyperplane arrangements. 
---
--the hyperplanes are in {z,t,w1,w2}

--fixHyperplaneRestriction=(z,t,wList)->SetHyperplaneIntersection=apply(#wList,i=>wList_i=>makeB'Section({z,t}|{1},NameB'Section=>wList_i));

















{*
--Example 1. 
restart
installPackage"BraidGroupHomotopy"
--
printingPrecision=300
R=CC[z,t,x,y,yA,s]
varList=(z,t,x,y,yA,s)
f=z^3-t; 
computeBranchPoints(z,t,f)
SetBraidGroupBranchPoints
cbpaf=computeBasePointAndFiber(z,t,f)
SetDownstairsStartPoint==first cbpaf
SetUpstairsStartFiber==last cbpaf
#SetDownstairsStartPoint==1
3==#SetUpstairsStartFiber
SetEncirclingTriangles=computeEncirclingTriangles(null)
#computeEncirclingTriangles(null)==#SetBraidGroupBranchPoints
oneTriangle=first SetEncirclingTriangles
--computeTwistLocusOfTriangle(varList,f,oneTriangle)
computeBraid(varList,f,oneTriangle)


i16 : --computeTwistLocusOfTriangle(varList,f,oneTriangle)
      computeBraid(varList,f,oneTriangle)
(Twist locus on a segment, (.277383318074136-.868247746428931*ii, -.0158990300749976-.047404860960371*ii), {.0811546636435964-.319040145184304*ii, .000741764701597314-.0939793584953007*ii, -1.11022302462516e-16-.0919032967007841*ii, -.0158990300749976-.047404860960371*ii})
(Twist locus on a segment, (-.0158990300749976-.047404860960371*ii, .0490033288920267+.0099334665397883*ii), {-.00989488052275946-.0421004637263497*ii, -8.67361737988404e-17-.0333587799259625*ii, .00133495180234849-.032179409794277*ii, .0177751671196669-.0176552158353456*ii, .0376787242322909-.0000713145055991904*ii, .0377594464967743+9.71445146547012e-17*ii, .0490033288920267+.0099334665397883*ii})
(Twist locus on a segment, (.0490033288920267+.0099334665397883*ii, -.0331042988171353+.0374713944206884*ii), {.0468752518895219+.010647198414098*ii, .0443049269897508+.0115092550061005*ii, .0409638259134696+.0126298208176785*ii, .0201452246906282+.0196121336068758*ii, 3.46944695195361e-17+.0263686038662816*ii, -.0000645877871833511+.0263902658468775*ii, -.0331042988171353+.0374713944206884*ii})
(Twist locus on a segment, (-.0331042988171353+.0374713944206884*ii, -.0158990300749976-.047404860960371*ii), {-.0319077651335191+.0315687084451196*ii, -.0315666884062127+.0298861241250145*ii, -.0297997877856911+.0211697297144156*ii, -.025538248583783+.000146896822753227*ii, -.0255084711235581, -.021223713731719-.0211373718180873*ii, -.0158990300749976-.047404860960371*ii})
(Twist locus on a segment, (-.0158990300749976-.047404860960371*ii, .277383318074136-.868247746428931*ii), {-2.42861286636753e-17-.0919032967007865*ii, .0292922168893675-.173886777756057*ii, .277383318074136-.868247746428931*ii})
Twist loci found! 
(END SEGMENT 0, .0811546636435964-.319040145184304*ii)
(END SEGMENT 1, .000741764701597314-.0939793584953007*ii)
The twist is : [3,2]
1 twists on fiber over: 
-.091903*ii
{-.20365, .10183+.17637*ii, .10183-.17637*ii}
{-.20365, .10183-.17637*ii, .10183+.17637*ii}
(END SEGMENT 2, -1.11022302462516e-16-.0919032967007841*ii)
(END SEGMENT 3, -.0158990300749976-.047404860960371*ii)
(END SEGMENT 4, -.00989488052275946-.0421004637263497*ii)
The twist is : [2,3]
1 twists on fiber over: 
-.033359*ii
{-.10363, .051814-.089744*ii, .051814+.089744*ii}
{-.10363, .051814+.089744*ii, .051814-.089744*ii}
(END SEGMENT 5, -8.67361737988404e-17-.0333587799259625*ii)
(END SEGMENT 6, .00133495180234849-.032179409794277*ii)
(END SEGMENT 7, .0177751671196669-.0176552158353456*ii)
(END SEGMENT 8, .0376787242322909-.0000713145055991904*ii)
The twist is : [1,2]
1 twists on fiber over: 
.037759
{-.056276-.097472*ii, -.056276+.097472*ii, .11255}
{-.056276+.097472*ii, -.056276-.097472*ii, .11255}
(END SEGMENT 9, .0377594464967743+9.71445146547012e-17*ii)
(END SEGMENT 10, .0490033288920267+.0099334665397883*ii)
(END SEGMENT 11, .0468752518895219+.010647198414098*ii)
(END SEGMENT 12, .0443049269897508+.0115092550061005*ii)
(END SEGMENT 13, .0409638259134696+.0126298208176785*ii)
(END SEGMENT 14, .0201452246906282+.0196121336068758*ii)
The twist is : [2,3]
1 twists on fiber over: 
.026369*ii
{-.088591, .044296-.076722*ii, .044296+.076722*ii}
{-.088591, .044296+.076722*ii, .044296-.076722*ii}
(END SEGMENT 15, 3.46944695195361e-17+.0263686038662816*ii)
(END SEGMENT 16, -.0000645877871833511+.0263902658468775*ii)
(END SEGMENT 17, -.0331042988171353+.0374713944206884*ii)
(END SEGMENT 18, -.0319077651335191+.0315687084451196*ii)
(END SEGMENT 19, -.0315666884062127+.0298861241250145*ii)
(END SEGMENT 20, -.0297997877856911+.0211697297144156*ii)
(END SEGMENT 21, -.025538248583783+.000146896822753227*ii)
The twist is : [1,2]
1 twists on fiber over: 
-.025508
{-.043327-.075045*ii, -.043327+.075045*ii, .086654}
{-.043327+.075045*ii, -.043327-.075045*ii, .086654}
(END SEGMENT 22, -.0255084711235581)
(END SEGMENT 23, -.021223713731719-.0211373718180873*ii)
(END SEGMENT 24, -.0158990300749976-.047404860960371*ii)
The twist is : [2,3]
1 twists on fiber over: 
-.091903*ii
{-.20365, .10183-.17637*ii, .10183+.17637*ii}
{-.20365, .10183+.17637*ii, .10183-.17637*ii}
(END SEGMENT 25, -2.42861286636753e-17-.0919032967007865*ii)
(END SEGMENT 26, .0292922168893675-.173886777756057*ii)
(END SEGMENT 27, .277383318074136-.868247746428931*ii)
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

--Example 3.    CHECK
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
computeTwistLocusOfTriangle(varList,f,oneTriangle)
computeBraid(varList,f,oneTriangle)


Twist loci found! 
(END SEGMENT 0, .707384468716571-.151197845733397*ii)
(END SEGMENT 1, .825186156838347-.107011869220404*ii)
The twist is : [2,3]
1 twists on fiber over: 
.96901-.05306*ii
{-.67295+1.1172*ii, -.63104-1.1414*ii, -.63104-1.0371*ii, -.58266+1.0651*ii, 1.2137-.02793*ii, 1.304+.0242*ii}
{-.67295+1.1172*ii, -.63104-1.0371*ii, -.63104-1.1414*ii, -.58266+1.0651*ii, 1.2137-.02793*ii, 1.304+.0242*ii}
(END SEGMENT 2, .969012506114898-.0530643597163281*ii)
(END SEGMENT 3, .984100969925038-.0474048609604068*ii)
(END SEGMENT 4, 1.02681016878938-.00967319649198254*ii)
(END SEGMENT 5, 1.03136629577763-.00564806234435318*ii)
The twist is : [2,1]
The twist is : [3,4]
The twist is : [5,6]
3 twists on fiber over: 
1.0378
{-.66589+1.0719*ii, -.66589-1.0719*ii, -.59534-1.1126*ii, -.59534+1.1126*ii, 1.2612-.04073*ii, 1.2612+.04073*ii}
{-.66589-1.0719*ii, -.66589+1.0719*ii, -.59534+1.1126*ii, -.59534-1.1126*ii, 1.2612+.04073*ii, 1.2612-.04073*ii}
(END SEGMENT 6, 1.03775944649685+2.77555756156289e-17*ii)
(END SEGMENT 7, 1.04900332889206+.00993346653975245*ii)
The twist is : [2,3]
1 twists on fiber over: 
.98096+.03276*ii
{-.66406-1.1116*ii, -.63062+1.0491*ii, -.63062+1.1309*ii, -.5932-1.0707*ii, 1.2238+.02161*ii, 1.2947-.01931*ii}
{-.66406-1.1116*ii, -.63062+1.1309*ii, -.63062+1.0491*ii, -.5932-1.0707*ii, 1.2238+.02161*ii, 1.2947-.01931*ii}
(END SEGMENT 8, .980955092695135+.0327560406792847*ii)
(END SEGMENT 9, .9668957011829+.0374713944206526*ii)
The twist is : [1,2]
The twist is : [4,3]
2 twists on fiber over: 
.97449
{-.6463-1.1194*ii, -.6463+1.1194*ii, -.61272+1.0613*ii, -.61272-1.0613*ii, 1.2254, 1.2926}
{-.6463+1.1194*ii, -.6463-1.1194*ii, -.61272-1.0613*ii, -.61272+1.0613*ii, 1.2254, 1.2926}
(END SEGMENT 10, .974491528876469-6.45317133063372e-16*ii)
(END SEGMENT 11, .978564023410521-.0200902462664542*ii)
(END SEGMENT 12, .978984039152504-.0221622489718406*ii)
The twist is : [2,3]
1 twists on fiber over: 
.98109-.03253*ii
{-.66395+1.1115*ii, -.63062-1.1308*ii, -.63062-1.0492*ii, -.59333+1.0707*ii, 1.224-.02153*ii, 1.2946+.01924*ii}
{-.66395+1.1115*ii, -.63062-1.0492*ii, -.63062-1.1308*ii, -.59333+1.0707*ii, 1.224-.02153*ii, 1.2946+.01924*ii}
(END SEGMENT 13, .981086135917627-.0325322178776618*ii)
(END SEGMENT 14, .984100969925038-.0474048609604068*ii)
The twist is : [3,2]
1 twists on fiber over: 
.96901-.05306*ii
{-.67295+1.1172*ii, -.63104-1.0371*ii, -.63104-1.1414*ii, -.58266+1.0651*ii, 1.2137-.02793*ii, 1.304+.0242*ii}
{-.67295+1.1172*ii, -.63104-1.1414*ii, -.63104-1.0371*ii, -.58266+1.0651*ii, 1.2137-.02793*ii, 1.304+.0242*ii}
(END SEGMENT 15, .969012506114882-.0530643597163341*ii)
(END SEGMENT 16, .959243929149076-.0567284338083932*ii)
(END SEGMENT 17, .929761676445676-.0677868672265756*ii)
(END SEGMENT 18, -.185372328453442-.486060032255219*ii)
(END SEGMENT 19, -.458553083967123-.588526801932193*ii)
(END SEGMENT 20, -.616600397378907-.647808421644363*ii)
END BRANCH POINT 
1234567
4126376

(1235764)
()
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
installPackage"Bertini"
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

SetDownstairsStartPoint={.05-0.05*ii}
SetUpstairsStartFiber=bertiniZeroDimSolve({sub(f,{t=>first SetDownstairsStartPoint})},AffVariableGroup=>{z}    )/coordinates//flatten//sort

#SetDownstairsStartPoint==1
#SetUpstairsStartFiber==first degree f
SetEncirclingTriangles=computeEncirclingTriangles(null)
#computeEncirclingTriangles(null)==#SetBraidGroupBranchPoints

theBraids=apply(SetEncirclingTriangles,oneTriangle->computeBraid(varList,f,oneTriangle))




----Example 5 \label{ex:dessinB}
restart
loadPackage"BraidGroupHomotopy"
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



{*

computeTwistLocusOfTriangle=method(TypicalValue=>Thing,Options=>{    	
    })
computeTwistLocusOfTriangle(Thing,Thing,Thing) := o ->(varList,f,aTriangle)->(
(z,t,x,y,yA,s):=varList;
oneTwistLocusPathAndLoop:={};
thePath:=SetDownstairsStartPoint|aTriangle|SetDownstairsStartPoint;
for segmentNumber to #thePath-1-1 do(
--print segmentNumber;
oneSegment:={thePath_segmentNumber,thePath_(segmentNumber+1)};
gam':=oneSegment_0;
gam'':=oneSegment_1;
fs:=sub(f,{z=>x+ii*y,	t=>(1-s)*gam'+s*gam''	});
--print 3;
g1:=value replace("ii","0",     toString  (fs));
--g1=1/sub(max((flatten entries ((coefficients g1)_1))),CC)*g1;
--print 1;
g2:=value replace("ii","0",     toString  (ii*fs));
--g2=1/sub(max((flatten entries ((coefficients g2)_1))),CC)*g2;
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
branchS:=delete(null,apply(realSCoords,i->if  (realPart i)<1 and 0<(realPart  i) then i else null));
if #branchS>0 then branchS=sort branchS;
twistLocusSegment:=    (for i in branchS list 	sub(sub((1-s)*gam'+s*gam'',{s=>i}),CC))|{gam''}    ;
print ("Twist locus on a segment",(gam',gam''),twistLocusSegment);
print("tw"|segmentNumber+1|" = "|toString twistLocusSegment);
oneTwistLocusPathAndLoop=oneTwistLocusPathAndLoop|twistLocusSegment;
--print (#oneTwistLocusPathAndLoop,segmentNumber)
);
print ("Twist loci found! ");
return oneTwistLocusPathAndLoop
)
*}