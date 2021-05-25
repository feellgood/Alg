a=1;
b=0.5;
H=0.5;

lc=0.005; // taille des mailles

Point(1) = { 0,  0, 0, lc};  // x, y, z=0, taille
Point(2) = { a,  0, 0, lc};
Point(3) = { 0,  b, 0, lc};
Point(4) = {-a,  0, 0, lc};
Point(5) = { 0, -b, 0, lc};

Point(6) = {-H,  0, 0, lc}; // coord de l'armature interne
Point(7) = {+H,  0, 0, lc};

Ellipse(1) = {2, 1, 3, 3}; // armature externe
Ellipse(2) = {3, 1, 4, 4};
Ellipse(3) = {4, 1, 5, 5};
Ellipse(4) = {5, 1, 2, 2};

Line(5) = {6, 1};
Line(6) = {1, 7};
Line(7) = {4, 6};
Line(8) = {7, 2};

Line Loop(9) = {1, 2, 7, 5, 6, 8}; // surface geometrique
Plane Surface(10) = {9}; 

Line Loop(11) = {7, 5, 6, 8, -4, -3}; // surface geometrique
Plane Surface(12) = {11}; 

Physical Surface(1) = {10, 12};  // region physique surfacique 1 = air
Physical Line(2) = {2, 1, 4, 3}; // region physique lineique 2   = armature ext
Physical Line(3) = {5, 6};       // region physique lineique 3   = armature int

