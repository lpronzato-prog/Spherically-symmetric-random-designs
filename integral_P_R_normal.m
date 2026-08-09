function PT = integral_P_R_normal(tt,sig,rr,d)
% function PT = integral_P_R_normal(tt,sig,rr,d)
%__________________________________________________________________________
% computes int I_v((d-1)/2,(d-1)/2)*phi(R) dR
% where R ~N(0,sig^2*eye(d)/d) and v=min(max((T.^2-(-r+R).^2)./(4*r*R),0),1);
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

sigma=@(t,R,r) min(max((t.^2-(-r+R).^2)./(4*r*R),0),1); % first parameter for beta inverse (truncated to [0,1])
psinormal=@(r,d) 1/(2^(d/2-1))*d^(d/2)/gamma(d/2)*r.^(d-1).*exp(-d*r.^2/2);
phinormal_design=@(r,d,sig) (1/sig)*psinormal(r/sig,d);
integral_P_R=@(t,sig,r,d) integral(@(R) betainc(sigma(t,R,r),(d-1)/2,(d-1)/2).*phinormal_design(R,d,sig),0,Inf);

PT=NaN*tt;
LT=length(tt); 
for i=1:LT
    t=tt(i);
    r=rr(i);
    PT(i)=integral_P_R(t,sig,r,d);
end    

end