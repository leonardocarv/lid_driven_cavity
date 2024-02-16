function [diffx,diffy] = doDiffusion(prx0,prx1,prx2,pry1,pry2,dx,dy)
diffy = (pry1 - 2*prx0 + pry2)/(dy*dy);
diffx = (prx1 - 2*prx0 + prx2)/(dx*dx);
