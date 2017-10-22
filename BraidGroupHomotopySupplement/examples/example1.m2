restart
path=prepend("/Users/jo/Documents/GitStuff/SymmetricMonodomy",path)
installPackage"SymmetricMonodromy"
--peek SymmetricMonodromy


--
  path=prepend("/Users/jo/Documents/GitStuff/SymmetricMonodomy",path)
installPackage"SymmetricMonodromy"
printingPrecision=300

  R=QQ[x,y,z]**QQ[u1,u2,u3];
  f1=x^5-u1;
  f2=y^6-u2;
  f3=z^3+u3
  fList={f1,f2,f3};
  xList={{x,y,z}};
  uList={{u1,u2,u3}}; 
  writeBertiniParameters(theDir,"start_parameters",{{1,1,-1}})
  writeStartFiberFile(theDir,{{1,1,1}})
  writeParameterHomotopyInputFile(theDir,xList,uList,fList)
  trackOverLoops(0,4,uList,"/Users/jo/Desktop/Dump/")
  unsortedFiber=importLoopOutput(0,4,xList,"/Users/jo/Desktop/Dump/")--An error catcher for incorrect variables would be useful.       
  S=CC[x,y,z]
  cleanFiber=cleanLoopOutput({gens S},x,{},unsortedFiber)
   sortedFiber=for i in cleanFiber/last list unsortedFiber_(i_0)_(i_1)
writeStartFiberFile(theDir,sortedFiber)
  trackOverLoops(6,10,uList,"/Users/jo/Desktop/Dump/")
  unsortedFiber=importLoopOutput(6,10,xList,"/Users/jo/Desktop/Dump/")--An error catcher for incorrect variables would be useful.       
 unsortedFiber/length
