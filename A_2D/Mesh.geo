SetFactory("OpenCASCADE");


x1 = 0;
x2 = 2*Pi;
n =1000;
lc1 = 0.08;
lc2 = 0.08;
lc3 = 0.08;
lc4 = 0.08;

//First interface Points and Lines
For i In {0:n-1}

  x = x1 + (x2 - x1) * i/n;

  r = 0.7512+0.0124*Cos(x)+0.0414*Cos(2*x)+0.0120*Cos(3*x)+0.0111*Cos(4*x)-0.0061*Cos(5*x)-0.0008*Cos(6*x)-0.00001*Cos(7*x)+0.0033*Cos(8*x)+ 0.0217*Sin(x)-0.0197*Sin(2*x)-0.0027*Sin(3*x)-0.0084*Sin(4*x)-0.0135*Sin(5*x)-0.0034*Sin(6*x)-0.0001*Sin(7*x)-0.0048*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc1}; //if you want the inner boundary mesh to be more refined change lc1 to a smaller number

EndFor

//Spline creates a continuous boundary using the points above
s1 = newreg;
Spline(s1) = {pList1[{0:n-1}],1};

//Second interface Points and Lines
For i In {n:2*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = 1.0284+0.0362*Cos(x)-0.0557*Cos(2*x)+0.0190*Cos(3*x)-0.0185*Cos(4*x)+0.0038*Cos(5*x)+0.0045*Cos(6*x)-0.0002*Cos(7*x)+0.0006*Cos(8*x)-0.1984*Sin(x)-0.0173*Sin(2*x)-0.0008*Sin(3*x)-0.0133*Sin(4*x)-0.0026*Sin(5*x)-0.0066*Sin(6*x)-0.0029*Sin(7*x)+0.0064*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc2};

EndFor

//Spline creates a continuous boundary using the points above
s2 = newreg;
Spline(s2) = {pList1[{n:2*n-1}],n+1};


//Third Interface Points and Lines
For i In {2*n:3*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = 1.2584+0.0362*Cos(x)-0.0557*Cos(2*x)+0.0190*Cos(3*x)-0.0185*Cos(4*x)+0.0038*Cos(5*x)+0.0045*Cos(6*x)-0.0002*Cos(7*x)+0.0006*Cos(8*x)-0.1984*Sin(x)-0.0173*Sin(2*x)-0.0008*Sin(3*x)-0.0133*Sin(4*x)-0.0026*Sin(5*x)-0.0066*Sin(6*x)-0.0029*Sin(7*x)+0.0064*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc3};

EndFor

//Spline creates a continuous boundary using the points above
s3 = newreg;
Spline(s3) = {pList1[{2*n:3*n-1}],2*n+1};



//Fourth interface Points and Lines
For i In {3*n:4*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = 1.3584+0.0362*Cos(x)-0.0557*Cos(2*x)+0.0190*Cos(3*x)-0.0185*Cos(4*x)+0.0038*Cos(5*x)+0.0045*Cos(6*x)-0.0002*Cos(7*x)+0.0006*Cos(8*x)-0.1984*Sin(x)-0.0173*Sin(2*x)-0.0008*Sin(3*x)-0.0133*Sin(4*x)-0.0026*Sin(5*x)-0.0066*Sin(6*x)-0.0029*Sin(7*x)+0.0064*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc4};

EndFor

//Spline creates a continuous boundary using the points above
s4 = newreg;
Spline(s4) = {pList1[{3*n:4*n-1}],3*n+1};

//The following lines should be written by hand. For extracting the labels you have to check them on GMSH. For example the number {2} for the first curveloop was 
//exracted by hovering the mouse pointer over the inner interface.

//+
Curve Loop(1)={2};
//+
Curve Loop(2)={1};
//+
Plane Surface(1) = {1,2};
//+
Curve Loop(3)={3};
//+
Curve Loop(4)={2};
//+
Plane Surface(2) = {3,4};
//+
Curve Loop(5)={4};
//+
Curve Loop(6)={3};
//+
Plane Surface(3) = {5,6};

//Surface and volume labeling
//+
Physical Curve(1) = {1};
//+
Physical Curve(2) = {2};
//+
Physical Curve(3) = {3};
//+
Physical Curve(4) = {4};
//+
Physical Surface(5) = {1};
//+
Physical Surface(6) = {2};
//+
Physical Surface(7) = {3};
