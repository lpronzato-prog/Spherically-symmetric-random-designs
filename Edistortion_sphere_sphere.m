function Edist = Edistortion_sphere_sphere(R,d,n,s)
% function Edist = Edistortion_sphere_sphere(R,d,n,s)
%__________________________________________________________________________
% Computes the expected distortion ( = (Ls-mean quantisation error)^s) for 
% a random design with n points independently uniformly distributed on the 
% sphere S_{d-1}(0,R) when the reference measure mu is uniform on the 
% sphere S_{d-1}(0,1).
% -- No approximation
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

r=1;
sigma=@(t,R) (t.^2-(-r+R).^2)./(4*r*R); % first parameter for beta inverse
integrant_Q=@(t,R,d,n,s) t.^(s-1).*( 1-betainc(sigma(t,R),(d-1)/2,(d-1)/2) ).^n;
% E{quantisation error^s}
Edist = s*integral(@(t) integrant_Q(t,R,d,n,s),1-R,1+R)+(1-R)^s; 

end