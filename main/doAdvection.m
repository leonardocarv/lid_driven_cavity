function [advx,advy] = doAdvection(psiu1,psiu2,psiv1,psiv2,prx1,prx2,pry1,pry2,dx,dy)
% this function compute the components of advection with a central finite
% difference approximation
advx = ((psiu1 - psiu2)/(2*dy))*((prx1 - prx2)/(2*dx));
advy = ((psiv1 - psiv2)/(2*dx))*((pry1 - pry2)/(2*dy));


