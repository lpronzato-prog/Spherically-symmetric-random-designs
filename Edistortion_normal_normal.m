function Edist = Edistortion_normal_normal(sig,d,n,s)
% function Edist = Edistortion_ball_normal(sig,d,n,s)
%__________________________________________________________________________
% Computes the expected distortion ( = (Ls-mean quantisation error)^s) for 
% a random design with n points independent normal N(0,sig^2*eye(d)/d)
% when the reference measure mu is normal N(0,eye(d)/d).
% -- No approximation
%
% Uses integral_P_R_normal.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

psinormal=@(r,d) 1/(2^(d/2-1))*d^(d/2)/gamma(d/2)*r.^(d-1).*exp(-d*r.^2/2);

% integral_P_R_normal(tt,sig,rr,d) computes 
%   integral_P_R=@(t,sig,r,d) integral(@(R) betainc(sigma(t,R,r),(d-1)/2,(d-1)/2).*phinormal_design(R,d,sig),0,Inf);
%   for all t=T(i) and r=rr(i), i=1,...,length(T)
% with
% sigma=@(t,R,r) min(max((t.^2-(-r+R).^2)./(4*r*R),0),1); % first parameter for beta inverse (truncated to [0,1])
% phinormal_design=@(r,d,sig) (1/sig)*psinormal(r/sig,d);

% E{quantisation error^s}
Edist =  s*integral2(@(t,r) t.^(s-1).*(1-integral_P_R_normal(t,sig,r,d)).^n.*psinormal(r,d),0,Inf,0,Inf) ; 

end