algebraic3d

solid cyl = cylinder(0,0,0; 0,0,5.000; 10.000) -maxh=0.500;
solid bigcyl = plane(0,0,0; 0,0,-1)
	 and plane(0,0,5.000; 0,0,1)
	 and cyl -bc=2;

tlo bigcyl; 
