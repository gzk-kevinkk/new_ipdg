algebraic3d

# example with two sub-domains

solid cube = plane (-1, -1, -1; 0, 0, -1)
         and plane (-1, -1, -1; 0, -1, 0)
         and plane (-1, -1, -1; -1, 0, 0)
         and plane (1, 1, 1; 0, 0, 1)
         and plane (1, 1, 1; 0, 1, 0)
         and plane (1, 1, 1; 1, 0, 0) -bc=1;
solid cutplane = plane (0.0, 0, 0; -1, 0, 0) -bc=4;

solid right = cube and cutplane;
solid left = cube and not cutplane;

tlo right -col=[1,0,0];
tlo left -col=[0,0,1];
