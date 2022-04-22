algebraic3d

# example with two sub-domains

solid cube = plane (0.0, 0, 0; 0, 0, -1)
         and plane (0.0, 0, 0; 0, -1, 0)
         and plane (0.0, 0, 0; -1, 0, 0)
         and plane (1, 1, 1; 0, 0, 1)
         and plane (1, 1, 1; 0, 1, 0)
         and plane (1, 1, 1; 1, 0, 0) -bc=1;

tlo cube -col=[1,0,0];
