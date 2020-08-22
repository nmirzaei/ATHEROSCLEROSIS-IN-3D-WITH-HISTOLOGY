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

  r = (0.714858)+(0.060788)*Cos(x)+(0.014752)*Cos(2*x)+(-0.006091)*Cos(3*x)+(0.004795)*Cos(4*x)+(0.000545)*Cos(5*x)+(-0.001939)*Cos(6*x)+(0.000596)*Cos(7*x)+(-0.002890)*Cos(8*x)+(0.079346)*Sin(x)+(-0.015798)*Sin(2*x)+(-0.007070)*Sin(3*x)+(-0.013298)*Sin(4*x)+(0.001152)*Sin(5*x)+(0.003392)*Sin(6*x)+(-0.001299)*Sin(7*x)+(0.001143)*Sin(8*x);


  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc1}; //if you want the inner boundary mesh to be more refined change lc1 to a smaller number

EndFor

s1 = newreg;
Spline(s1) = {pList1[{0:n-1}],1};

//Second interface Points and Lines

For i In {n:2*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = (1.034536)+(-0.030658)*Cos(x)+(-0.066222)*Cos(2*x)+(0.010444)*Cos(3*x)+(-0.004373)*Cos(4*x)+(0.002543)*Cos(5*x)+(-0.000342)*Cos(6*x)+(0.001576)*Cos(7*x)+(-0.000734)*Cos(8*x)+(-0.163745)*Sin(x)+(0.042981)*Sin(2*x)+(-0.006253)*Sin(3*x)+(0.000991)*Sin(4*x)+(0.005113)*Sin(5*x)+(-0.001116)*Sin(6*x)+(-0.004498)*Sin(7*x)+(0.002001)*Sin(8*x);


  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc2};

EndFor

s2 = newreg;
Spline(s2) = {pList1[{n:2*n-1}],n+1};


//Third Interface Points and Lines

For i In {2*n:3*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = (1.234536)+(-0.030658)*Cos(x)+(-0.066222)*Cos(2*x)+(0.010444)*Cos(3*x)+(-0.004373)*Cos(4*x)+(0.002543)*Cos(5*x)+(-0.000342)*Cos(6*x)+(0.001576)*Cos(7*x)+(-0.000734)*Cos(8*x)+(-0.163745)*Sin(x)+(0.042981)*Sin(2*x)+(-0.006253)*Sin(3*x)+(0.000991)*Sin(4*x)+(0.005113)*Sin(5*x)+(-0.001116)*Sin(6*x)+(-0.004498)*Sin(7*x)+(0.002001)*Sin(8*x);



  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc3};

EndFor


s3 = newreg;
Spline(s3) = {pList1[{2*n:3*n-1}],2*n+1};



//Fourth interface Points and Lines

For i In {3*n:4*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = (1.334536)+(-0.030658)*Cos(x)+(-0.066222)*Cos(2*x)+(0.010444)*Cos(3*x)+(-0.004373)*Cos(4*x)+(0.002543)*Cos(5*x)+(-0.000342)*Cos(6*x)+(0.001576)*Cos(7*x)+(-0.000734)*Cos(8*x)+(-0.163745)*Sin(x)+(0.042981)*Sin(2*x)+(-0.006253)*Sin(3*x)+(0.000991)*Sin(4*x)+(0.005113)*Sin(5*x)+(-0.001116)*Sin(6*x)+(-0.004498)*Sin(7*x)+(0.002001)*Sin(8*x);


  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc4};

EndFor


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
