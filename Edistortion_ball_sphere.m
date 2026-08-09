function Edist = Edistortion_ball_sphere(R,d,n,s)
% function Edist = Edistortion_ball_sphere(R,d,n,s)
%__________________________________________________________________________
% Computes the expected distortion ( = (Ls-mean quantisation error)^s) for 
% a random design with n points independently uniformly distributed on the 
% sphere S_{d-1}(0,R) when the reference measure mu is uniform in B_d(0,1).
% -- No approximation
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

sigma=@(t,R,r) min(max((t.^2-(-r+R).^2)./(4*r*R),0),1); % first parameter for beta inverse (truncated to [0,1])
one_minuP_to_n=@(t,R,r,d,n) ( 1-betainc(sigma(t,R,r),(d-1)/2,(d-1)/2) ).^n;
psiball=@(r,d) d*r.^(d-1); % density of ||U|| for U uniform in B_d(0,1)

% E{quantisation error^s}
Edist = s*integral2(@(t,r) t.^(s-1).*one_minuP_to_n(t,R,r,d,n).*psiball(r,d),0,Inf,0,1);

end