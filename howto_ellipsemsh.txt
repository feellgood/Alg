1)  gmsh ellipse.geo -2 -format msh22 -o ellipseHRes.msh
2) rajouter les paramètres à la fin du .msh ( + une ligne vide après labalise $EndMeshFormat)

$Scale 
1.0

$Parameters 
epsr    1       1.  
rho     1       0.
$EndParameters

$Conditions 
V	2   0.  
V	3   1.
$EndDirichlet

3) ./exe ellipseHRes.msh
