SetFactory("OpenCASCADE");


x1 = 0;
x2 = 2*Pi;
n =1000;
lc1 = 0.04;
lc2 = 0.04;
lc3 = 0.04;
lc4 = 0.04;

//First interface Points and Lines

For i In {0:n-1}

  x = x1 + (x2 - x1) * i/n;

  r = (0.619235)+(-0.308003)*Cos(x)+(0.001217)*Cos(2*x)+(0.021283)*Cos(3*x)+(-0.012021)*Cos(4*x)+(-0.000073)*Cos(5*x)+(0.001816)*Cos(6*x)+(0.001903)*Cos(7*x)+(0.001001)*Cos(8*x)+(0.086732)*Sin(x)+(-0.063993)*Sin(2*x)+(0.008916)*Sin(3*x)+(-0.000986)*Sin(4*x)+(0.004452)*Sin(5*x)+(-0.002024)*Sin(6*x)+(0.000823)*Sin(7*x)+(-0.000632)*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc1}; //if you want the inner boundary mesh to be more refined change lc1 to a smaller number

EndFor

s1 = newreg;
Spline(s1) = {pList1[{0:n-1}],1};

//Second interface Points and Lines

For i In {n:2*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = (1.023553)+(-0.018985)*Cos(x)+(-0.020950)*Cos(2*x)+(0.013332)*Cos(3*x)+(-0.014868)*Cos(4*x)+(-0.002363)*Cos(5*x)+(0.000426)*Cos(6*x)+(0.003187)*Cos(7*x)+(-0.003068)*Cos(8*x)+(-0.045158)*Sin(x)+(-0.043990)*Sin(2*x)+(0.010668)*Sin(3*x)+(0.004061)*Sin(4*x)+(0.003894)*Sin(5*x)+(-0.001746)*Sin(6*x)+(0.001121)*Sin(7*x)+(-0.003120)*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc2};

EndFor

s2 = newreg;
Spline(s2) = {pList1[{n:2*n-1}],n+1};


//Third Interface Points and Lines

For i In {2*n:3*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = (1.223553)+(-0.018985)*Cos(x)+(-0.020950)*Cos(2*x)+(0.013332)*Cos(3*x)+(-0.014868)*Cos(4*x)+(-0.002363)*Cos(5*x)+(0.000426)*Cos(6*x)+(0.003187)*Cos(7*x)+(-0.003068)*Cos(8*x)+(-0.045158)*Sin(x)+(-0.043990)*Sin(2*x)+(0.010668)*Sin(3*x)+(0.004061)*Sin(4*x)+(0.003894)*Sin(5*x)+(-0.001746)*Sin(6*x)+(0.001121)*Sin(7*x)+(-0.003120)*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc3};

EndFor


s3 = newreg;
Spline(s3) = {pList1[{2*n:3*n-1}],2*n+1};



//Fourth interface Points and Lines

For i In {3*n:4*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = (1.323553)+(-0.018985)*Cos(x)+(-0.020950)*Cos(2*x)+(0.013332)*Cos(3*x)+(-0.014868)*Cos(4*x)+(-0.002363)*Cos(5*x)+(0.000426)*Cos(6*x)+(0.003187)*Cos(7*x)+(-0.003068)*Cos(8*x)+(-0.045158)*Sin(x)+(-0.043990)*Sin(2*x)+(0.010668)*Sin(3*x)+(0.004061)*Sin(4*x)+(0.003894)*Sin(5*x)+(-0.001746)*Sin(6*x)+(0.001121)*Sin(7*x)+(-0.003120)*Sin(8*x);


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
