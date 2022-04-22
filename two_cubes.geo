algebraic3d

# example with two sub-domains

solid cube1 = plane (0, 0, 0; 0, 0, -1)
         and plane (0, 0, 0; 0, -1, 0)
         and plane (0, 0, 0; -1, 0, 0)
         and plane (1, 1, 1; 0, 0, 1)
         and plane (1, 1, 1; 0, 1, 0)
         and plane (1, 1, 1; 1, 0, 0) -bc=1;
solid cube2 = plane (0.25, 0.25, 0.25; 0, 0, -1)
         and plane (0.25, 0.25, 0.25; 0, -1, 0)
         and plane (0.25, 0.25, 0.25; -1, 0, 0)
         and plane (0.75, 0.75, 0.75; 0, 0, 1)
         and plane (0.75, 0.75, 0.75; 0, 1, 0)
         and plane (0.75, 0.75, 0.75; 1, 0, 0) -bc=5;

solid in = cube2;
solid out = cube1 and not cube2;

tlo in -col=[1,0,0];
tlo out -col=[0,0,1];