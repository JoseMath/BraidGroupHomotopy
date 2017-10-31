restart
loadPackage"BraidGroupHomotopy"
--
printingPrecision=300
R=CC[z,t,x,y,yA,s]
varList=(z,t,x,y,yA,s)
f=z^4-4*z^2+3+t; 
RadiusBranchPoint=.5
SetBraidGroupBranchPoints={-3,1};
SetUpstairsStartFiber={-1.77273369086864+.0749433780110061*ii, -.939691646423518-.141381113276106*ii, .939691646423518+.141381113276106*ii, 1.77273369086864-.0749433780110058*ii}
SetDownstairsStartPoint={-.222095320092058+.604206069448226*ii}

#SetUpstairsStartFiber==first degree f
SetEncirclingTriangles=computeEncirclingTriangles(null)
#computeEncirclingTriangles(null)==#SetBraidGroupBranchPoints
oneTriangle=first SetEncirclingTriangles
twoTriangle=last SetEncirclingTriangles

--computeTwistLocusOfTriangle(varList,f,oneTriangle)
computeBraid(varList,f,oneTriangle)

matLabTriangle=aTri->(
    print ("x = "|toString new Array from  ((SetDownstairsStartPoint|aTri|SetDownstairsStartPoint)/realPart));
print ("y= "|toString new Array from  ((SetDownstairsStartPoint|aTri|SetDownstairsStartPoint)/imaginaryPart)
    ))

matLabTriangle(oneTriangle)
matLabTriangle(twoTriangle)

-3.2551   [2,3]

computeBraid(varList,f,twoTriangle)
1.3776 [1,2][3,4]
    
    end
    ---
    Matlab code to plot
    
x1 = [-.222095320092058, -2.50996671107938, -3.331042988171, -3.15899030074962, -2.50996671107938, -.222095320092058]
x2 = [-.222095320092058, .668957011829001, .841009699250378, 1.49003328892062, .668957011829001, -.222095320092058]
y1= [.604206069448226, .0993346653975306, .374713944206532, -.474048609604062, .0993346653975306, .604206069448226]
y2=[    .604206069448226, .374713944206532, -.474048609604062, .0993346653975306, .374713944206532, .604206069448226]

xB1=[1]
yB1=[0]
xB2=[-3]
yB2=[0]

Tx=[-3.2551,1.3776]
Ty=[0,0]
[2,3]
[1,2][3,4]


plot (x1, y1,'b',x2,y2,'b',xB1,yB1,'r*',xB2,yB2,'r*',Tx,Ty,'mo')
help plot

The twist is : [2,3]
1 twists on fiber over: 
-3.2551
{-2.0156, -.25057*ii, .25057*ii, 2.0156}
{-2.0156, .25057*ii, -.25057*ii, 2.0156}



The twist is : [1,2]
The twist is : [3,4]
2 twists on fiber over: 
{-1.4304-.21479*ii, -1.4304+.21479*ii, 1.4304-.21479*ii, 1.4304+.21479*ii}
{-1.4304+.21479*ii, -1.4304-.21479*ii, 1.4304+.21479*ii, 1.4304-.21479*ii}










