function Ropt = Ropt_sphere_sphere(d,n,s,accuracy)
% function Ropt = Ropt_sphere_sphere(d,n,s,accuracy)
%__________________________________________________________________________
% Computes the optimal radius Ropt that minimizes the expected distortion 
% ( = (Ls-mean quantisation error)^s) for a random design with n points 
% independently uniformly distributed on the sphere S_{d-1}(0,R) when the 
% reference measure mu is uniform on the sphere S_{d-1}(0,1)
% accuracy = precision required on Ropt (e.g., 1e-4)
%
% Uses Edistortion_sphere_sphere.m, GoldenS.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

% E{quantisation error^s} for n points uniformly distributed in the ball or on the sphere
ED_R=@(R) Edistortion_sphere_sphere(R,d,n,s);
[aN,bN]=GoldenS(ED_R,0.1,1,accuracy);  
Ropt=(aN+bN)/2;
end
