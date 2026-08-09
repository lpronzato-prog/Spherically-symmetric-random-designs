function Ropt = Ropt_normal(design,d,n,s,accuracy)
% function Ropt = Ropt_normal(design,d,n,s,accuracy)
%__________________________________________________________________________
% Computes the optimal std sigopt or radius Ropt that minimises 
% the expected distortion ( = (Ls-mean quantisation error)^s) 
% for a random design with n points independently  
%   - normally distributed N(0,sig^2*eye(d)/d) when design = 'normal'
%   - uniformly distributed on the sphere S_{d-1}(0,R) when design = 'sphere'
% for the reference measure mu normal N(0,eye(d)/d)
% accuracy = precision required on Ropt (e.g., 1e-4)
%
% Uses Edistortion_normal_normal.m, Edistortion_normal_sphere.m, GoldenS.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

% E{quantisation error^s} for n points uniformly distributed in the ball or on the sphere
switch lower(design)
    case 'normal'
        ED_R=@(R) Edistortion_normal_normal(R,d,n,s);
    case 'sphere'
        ED_R=@(R) Edistortion_normal_sphere(R,d,n,s);
    end    
[aN,bN]=GoldenS(ED_R,0.1,1,accuracy);  
Ropt=(aN+bN)/2;
end
