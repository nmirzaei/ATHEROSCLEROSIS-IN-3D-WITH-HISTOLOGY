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

  r = (0.578244)+(-0.372715)*Cos(x)+(0.034017)*Cos(2*x)+(0.011333)*Cos(3*x)+(0.006484)*Cos(4*x)+(-0.005619)*Cos(5*x)+(-0.004827)*Cos(6*x)+(0.004946)*Cos(7*x)+(0.000575)*Cos(8*x)+(0.003674)*Sin(x)+(-0.023116)*Sin(2*x)+(0.012197)*Sin(3*x)+(0.005207)*Sin(4*x)+(0.004428)*Sin(5*x)+(0.002318)*Sin(6*x)+(-0.005940)*Sin(7*x)+(0.004274)*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc1}; //if you want the inner boundary mesh to be more refined change lc1 to a smaller number

EndFor

s1 = newreg;
Spline(s1) = {pList1[{0:n-1}],1};

//Second interface Points and Lines

For i In {n:2*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = (1.042974)+(-0.058644)*Cos(x)+(0.001808)*Cos(2*x)+(-0.013822)*Cos(3*x)+(0.003230)*Cos(4*x)+(-0.006332)*Cos(5*x)+(-0.003010)*Cos(6*x)+(0.000307)*Cos(7*x)+(0.001640)*Cos(8*x)+(-0.031070)*Sin(x)+(-0.000964)*Sin(2*x)+(-0.003722)*Sin(3*x)+(-0.007314)*Sin(4*x)+(0.001668)*Sin(5*x)+(-0.001127)*Sin(6*x)+(0.005961)*Sin(7*x)+(0.004657)*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc2};

EndFor

s2 = newreg;
Spline(s2) = {pList1[{n:2*n-1}],n+1};


//Third Interface Points and Lines

For i In {2*n:3*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = (1.242974)+(-0.058644)*Cos(x)+(0.001808)*Cos(2*x)+(-0.013822)*Cos(3*x)+(0.003230)*Cos(4*x)+(-0.006332)*Cos(5*x)+(-0.003010)*Cos(6*x)+(0.000307)*Cos(7*x)+(0.001640)*Cos(8*x)+(-0.031070)*Sin(x)+(-0.000964)*Sin(2*x)+(-0.003722)*Sin(3*x)+(-0.007314)*Sin(4*x)+(0.001668)*Sin(5*x)+(-0.001127)*Sin(6*x)+(0.005961)*Sin(7*x)+(0.004657)*Sin(8*x);

  pList1[i] = newp;

  Point (pList1[i]) = {r*Cos (x),r*Sin (x), -1, lc3};

EndFor


s3 = newreg;
Spline(s3) = {pList1[{2*n:3*n-1}],2*n+1};



//Fourth interface Points and Lines

For i In {3*n:4*n-1}

  x = x1 + (x2 - x1) * i/n;

  r = (1.342974)+(-0.058644)*Cos(x)+(0.001808)*Cos(2*x)+(-0.013822)*Cos(3*x)+(0.003230)*Cos(4*x)+(-0.006332)*Cos(5*x)+(-0.003010)*Cos(6*x)+(0.000307)*Cos(7*x)+(0.001640)*Cos(8*x)+(-0.031070)*Sin(x)+(-0.000964)*Sin(2*x)+(-0.003722)*Sin(3*x)+(-0.007314)*Sin(4*x)+(0.001668)*Sin(5*x)+(-0.001127)*Sin(6*x)+(0.005961)*Sin(7*x)+(0.004657)*Sin(8*x);


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
