algebraic3d

#define box:
#solid box = orthobrick (0, , 0; 40, 10,10);
solid outerborder = cylinder(0,0,0; 0,0,4.0; 2.50)
		    and plane(0,0,0;0,0,-1)
		    and plane(0,0,4.0;0,0,1);



# define electrode bars

##############  bar1 #################

######## Monitoring Wells ###########

solid cyl1 = cylinder (1.25, 0, 0; 1.25, 0, 1.0; 0.05) 
	and plane (0, 0, 0; 0, 0, -1)
	and plane (0, 0, 1.0; 0, 0, 1);
	
solid cyl2 = cylinder (1.25, 0, 1.0; 1.25, 0, 1.1; 0.05)
	and plane (0, 0, 1.0; 0, 0, -1)
	and plane (0, 0, 1.1; 0, 0, 1);

solid cyl3 = cylinder (1.25, 0, 1.1; 1.25, 0, 1.271; 0.05)
	and plane(0, 0, 1.1; 0, 0, -1)
	and plane(0, 0, 1.271; 0, 0, 1);

solid cyl4 = cylinder (1.25, 0, 1.271; 1.25, 0, 1.371; 0.05) 
	and plane (0, 0, 1.271; 0, 0, -1)
	and plane (0, 0, 1.371; 0, 0, 1);
	
solid cyl5 = cylinder (1.25, 0, 1.371; 1.25, 0, 1.542; 0.05)
	and plane (0, 0, 1.371; 0, 0, -1)
	and plane (0, 0, 1.542; 0, 0, 1);

solid cyl6 = cylinder (1.25, 0, 1.542; 1.25, 0, 1.642; 0.05)
	and plane (0, 0, 1.542; 0, 0, -1)
	and plane (0, 0, 1.642; 0, 0, 1);

solid cyl7 = cylinder (1.25, 0, 1.642; 1.25, 0, 1.813; 0.05)
	and plane (0, 0, 1.642; 0, 0, -1)
	and plane (0, 0, 1.813; 0, 0, 1);

solid cyl8 = cylinder (1.25, 0, 1.813; 1.25, 0, 1.913; 0.05)
	and plane (0, 0, 1.813; 0, 0, -1)
	and plane (0, 0, 1.913; 0, 0, 1);

solid cyl9 = cylinder (1.25, 0, 1.913; 1.25, 0, 2.084; 0.05)
	and plane (0, 0, 1.913; 0, 0, -1)
	and plane (0, 0, 2.084; 0, 0, 1);

solid cyl10 = cylinder (1.25, 0, 2.084; 1.25, 0, 2.184; 0.05)
	and plane (0, 0, 2.084; 0, 0, -1)
	and plane (0, 0, 2.184; 0, 0, 1);

solid cyl11 = cylinder (1.25, 0, 2.184; 1.25, 0, 2.355; 0.05)
	and plane (0, 0, 2.184; 0, 0, -1)
	and plane (0, 0, 2.355; 0, 0, 1);

solid cyl12 = cylinder (1.25, 0, 2.355; 1.25, 0, 2.455; 0.05)
	and plane (0, 0, 2.355; 0, 0, -1)
	and plane (0, 0, 2.455; 0, 0, 1);

solid cyl13 = cylinder (1.25, 0, 2.455; 1.25, 0, 2.626; 0.05)
	and plane (0, 0, 2.455; 0, 0, -1)
	and plane (0, 0, 2.626; 0, 0, 1);

solid cyl14 = cylinder (1.25, 0, 2.626; 1.25, 0, 2.726; 0.05)
	and plane (0, 0, 2.626; 0, 0, -1)
	and plane (0, 0, 2.726; 0, 0, 1);

solid cyl15 = cylinder (1.25, 0, 2.726; 1.25, 0, 2.897; 0.05)
	and plane (0, 0, 2.726; 0, 0, -1)
	and plane (0, 0, 2.897; 0, 0, 1);

solid cyl16 = cylinder (1.25, 0, 2.897; 1.25, 0, 2.997; 0.05)
	and plane (0, 0, 2.897; 0, 0, -1)
	and plane (0, 0, 2.997; 0, 0, 1);

solid cyl17 = cylinder (1.25, 0, 2.997; 1.25, 0, 3.5; 0.05)
	and plane (0, 0, 2.997; 0, 0, -1)
	and plane (0, 0, 3.5; 0, 0, 1);


solid bar1 =cyl1 or cyl2 or cyl3 or cyl4 or cyl5 or cyl6 or
            cyl7 or cyl8 or cyl9 or cyl10 or cyl11 or
            cyl12 or cyl13 or cyl14 or cyl15 or cyl16 or
	    cyl17; 


#make copies

#solid cylm1 = multitranslate(-1.0,1.0,0;1; bar1);
#solid cylm = multitranslate(-1.0,-1.0,0;1; cylm1);


############ Geophysical Holes ###########

solid bar2 = cylinder (0, 0, 0; 0, 0, 3.5; 0.05)
        and plane (0, 0, 0; 0, 0, -1)
        and plane (0, 0, 3.5; 0, 0, 1);

solid cylg1 = bar1; #multitranslate(0.25,0,0;1; bar1);
solid cylg2 = multitranslate(-0.3661,0.8839,0;1; bar1);
solid cylg3 = multitranslate(-1.25,1.25,0;1; bar1);
solid cylg4 = multitranslate(-2.1339,0.8839,0;1; bar1);
solid cylg5 = multitranslate(-2.5,0,0;1; bar1);
solid cylg6 = multitranslate(-2.1339,-0.8839,0;1; bar1);
solid cylg7 = multitranslate(-1.25,-1.25,0;1; bar1);
solid cylg8 = multitranslate(-0.3661,-0.8839,0;1; bar1);

solid cylg  = cylg1 or cylg2 or cylg3 or cylg4 or cylg5 or
              cylg6 or cylg7 or cylg8;



############# Injection Well #############

#solid cyli = bar2; #multitranslate(0,-1.25,0;1; bar2);


######### combine ######

#solid main = outerborder and not cylm and not cylg
#	     and not cyli;

solid main = outerborder and not cylg;

#define sub-domains:
tlo main;

#boundingbox (-1, -1, -1; 11, 11, 11);

