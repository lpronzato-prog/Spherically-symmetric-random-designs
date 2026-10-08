function [Radii, Weights,ED,Sensitivities_on_support,Sensitivities_on_RHO] = Optimal_Sphere_Mixture(distribution,d,n,s,RHO,epsilon,k_VD,k_VE,k_remove,small_weight,k_merge,dist_min,initial)
% function [Radii, Weights,ED,Sensitivities_on_support,Sensitivities_on_RHO] = Optimal_Sphere_Mixture(distribution,d,n,s,RHO,epsilon,k_VD,k_VE,k_remove,small_weight,k_merge,dist_min,initial)
%  Computes the optimal distribution, in the form of a discrete measure
%  allocating the weights in the vector Weigths to the radii in Radii, for
%  the minimisation of the expected distortion of random designs for the
%  target measure mu
% -- uniform in the ball B_d(0,1) if distribution = 'ball'
% -- normal N(0,eye(d)/d) if distribution = 'normal'
% The algorithm uses both a vertex-direction and vertex-exchange method;
% see Section 9.1 of [Pronzato & Pazman, Springer, 2013]
% RHO (a row vector) specifies a set of candidate support points (e.g.,
% RH0=0:0.01:1 for the ball, RH0=0:0.03:2 for the normal distribution), 
% but optimal radii are also searched outside of this set (in its convex
% hull, however).
% The function returns the row vectors Radii and Weights that define the
% optimal distribution, together with 
% -- ED, the expected distortion for the optimal distribution,
% -- Sensitivities_on_support, the vector of sensitivities for R in Radii, and
% -- Sensitivities_on_RHO, the sensitivities computed for values of rho in RHO. 
% The algorithm stops when the efficiency (ratio of the optimal (unknown)
% distortion to the current one) is larger than 1-epsilon.
% k_VD specifies how often a vertex-direction step is used
% k_VE specifies how often a vertex-exchange step is used
% k_remove specifies how often points with weigth < small_weight are
%   potentially removed
% k_merge specifies how often points whose distance is less than dist_min
%   are potentially merged
% Typical values are k_VD=k_VE=k_merge=1, k_remove=5, small_weight=1e-3,
%   dist_min=0.2. The value of epsilon should be adapted to the difficulty of
%   the problem, values about 1e-4 or 1e-5 are often satisfactory. 
% When initial = '1pt-opt', the algorithm is initialised at the optimal
%   one-point measure, otherwise it is initalised at mean(RHO)
%-----
% Author: L. Pronzato, 2026 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

LRHO_0=length(RHO);
sigma=@(t,R,r) min(max((t.^2-(-r+R).^2)./(4*r*R),0),1); % first parameter for beta inverse (truncated to [0,1])
switch lower(distribution)
    case 'ball'
        psi=@(r,d) d*r.^(d-1);                  % U uniform in B_d(0,1)
        UB=1;                                   % upper bound for integration
        UB_GS=1;                                % upper bound for GS algorithm
    case 'normal'
        psi=@(r,d) 1/(2^(d/2-1))*d^(d/2)/gamma(d/2)*r.^(d-1).*exp(-d*r.^2/2); % U normal N(0,Id/d)
        UB=Inf;                                 % upper bound for integration
        UB_GS=2;                                % upper bound for GS algorithm
end
one_minuP_to_n=@(t,R,r,d,n) ( 1-betainc(sigma(t,R,r),(d-1)/2,(d-1)/2) ).^n;
% Expected distortion for n random points uniform on S_{d-1}(a);
ED=@(R,d,n,s) s*integral2(@(t,r) t.^(s-1).*one_minuP_to_n(t,R,r,d,n).*psi(r,d),0,Inf,0,UB); % E{distortion}
A_delta_R=@(t,R,r) betainc(sigma(t,R,r),(d-1)/2,(d-1)/2);

k=1;
Rk=mean(RHO);
switch lower(initial)
    case '1pt-opt'
    % initialisation by the optimal sphere
    ED_a=@(R) ED(R,d,n,s);
    [aN,bN]=GoldenS(ED_a,1e-2,UB_GS,1e-4); 
    a_star=(aN+bN)/2; 
    Rk=a_star;
end    
Rk
[RHO,~,ic]=unique([RHO Rk]);
LRHO=length(RHO);
Indices=[ic(end)]; % the last one, unless Rk was already in RHO
Radii=[Rk];
Weights=[1];
A_Phik=@(t,r) A_delta_R(t,Rk,r);
% E distortion
EDk=s*integral2(@(t,r) t.^(s-1).*(1-A_Phik(t,r)).^n .*psi(r,d),0,Inf,0,UB); % E{distortion}
% sensitivity
S_Phik=@(rho) n*s*integral2(@(t,r) t.^(s-1).*(1-A_Phik(t,r)).^(n-1) .*A_delta_R(t,rho,r) .*psi(r,d),0,Inf,0,UB);
S_bar_Phik=S_Phik(Rk)
% check optimality
Sens=NaN(1,LRHO);
for i=1:LRHO
    rho=RHO(i);
    Sens(i)=S_Phik(rho);
end
one_minusEff=(max(Sens-S_bar_Phik))/EDk;

while one_minusEff>epsilon   
    k=k+1;
    [k length(RHO)]
    EDkold=EDk;
    if rem(k,k_VD)==0 % vertex direction
        % optimal direction
        [~,imax]=max(Sens);
        Rk=RHO(imax); Ik=imax;
        A_Phikp1=@(alpha,t,r) (1-alpha)*A_Phik(t,r)+alpha*A_delta_R(t,Rk,r);
        ED_Phikp1=@(alpha) s*integral2(@(t,r) t.^(s-1).*(1-A_Phikp1(alpha,t,r)).^n .*psi(r,d),0,Inf,0,UB);
        [aN,bN]=GoldenS(ED_Phikp1,0,1,1e-6); 
        alpha_star=(aN+bN)/2;
        [lia,locb]=ismember(Ik,Indices);
        if lia==1
            % this radius is already there
            Weights=(1-alpha_star)*Weights;
            Weights(locb)=Weights(locb)+alpha_star;
        else
            % this is a new radius
            Indices=[Indices Ik];
            Radii=[Radii Rk];
            Weights=[(1-alpha_star)*Weights alpha_star];
        end 
        EDk=ED_Phikp1(alpha_star);
        A_Phik=@(t,r) Weights(1)*A_delta_R(t,Radii(1),r);
        if length(Indices)>1
            for i=2:length(Indices)
                A_Phik=@(t,r) A_Phik(t,r)+Weights(i)*A_delta_R(t,Radii(i),r);
            end
        end
    end
    if rem(k,k_VE)==0 % vertex exchange
        % optimal direction: Rkp-Rkm
        [~,imax]=max(Sens);
        [~,imink]=min(Sens(Indices)); 
        Rkp=RHO(imax); Ikp=imax;
        Rkm=Radii(imink); 
        A_Phikp1=@(alpha,t,r) A_Phik(t,r)+alpha*(A_delta_R(t,Rkp,r)-A_delta_R(t,Rkm,r));
        ED_Phikp1=@(alpha) s*integral2(@(t,r) t.^(s-1).*(1-A_Phikp1(alpha,t,r)).^n .*psi(r,d),0,Inf,0,UB);
        alpha_max=Weights(imink);
        [aN,bN]=GoldenS(ED_Phikp1,0,alpha_max,1e-6); 
        alpha_star=(aN+bN)/2;
        [lia,locb]=ismember(Ikp,Indices);
        if lia==1
            % this radius is already there
            Weights(locb)=Weights(locb)+alpha_star;
            Weights(imink)=Weights(imink)-alpha_star;
        else
            % this is a new radius
            Indices=[Indices Ikp];
            Radii=[Radii Rkp];
            Weights(imink)=Weights(imink)-alpha_star;
            Weights=[Weights alpha_star];
        end
        EDk=ED_Phikp1(alpha_star);
        A_Phik=@(t,r) Weights(1)*A_delta_R(t,Radii(1),r);
        if length(Indices)>1
            for i=2:length(Indices)
                A_Phik=@(t,r) A_Phik(t,r)+Weights(i)*A_delta_R(t,Radii(i),r);
            end
        end
    end
    
    if rem(k,k_remove)==0 % only try to remove points with negligible weights (only from time to time)
        % sort by increasing weights
        [~,isort]=sort(Weights,'ascend');
        Weights=Weights(isort);
        Radii=Radii(isort);
        Indices=Indices(isort);
        [~,i_small]=find(Weights<small_weight);
        while isempty(i_small)~=1 % try to remove small weights one by one
            i_remove=i_small(1);
            i_small(1)=[]; % removed or not, this one has already been considered           
            wmin=Weights(i_remove);
            A_Phik_try=@(t,r) (A_Phik(t,r)-wmin*A_delta_R(t,Radii(i_remove),r))/(1-wmin);
            EDk_try = s*integral2(@(t,r) t.^(s-1).*(1-A_Phik_try(t,r)).^n .*psi(r,d),0,Inf,0,UB); % E{distortion}
            if EDk_try<=(1+1e-6)*EDk % *(1+1e-6) to account for numerical inacurracy (we favour designs with small support)
                % remove the point
                A_Phik=@(t,r) A_Phik_try(t,r);
                EDk = EDk_try;
                Radii(i_remove)=[];
                Indices(i_remove)=[];
                Weights(i_remove)=[];
                Weights=Weights/sum(Weights);
            end
        end
        
    end

    if rem(k,k_merge)==0 % try merge radii (only from time to time)
        % sort the radii by increasing values (and reorder Weights and Indices accordingly)
        [Radii,isort]=sort(Radii,'ascend');
        Weights=Weights(isort);
        Indices=Indices(isort);
    
        % try to merge points close to each other
        [dmin,iclose]=min(diff(Radii)); % Radii(iclose) and Radii(iclose+1) are the closest ones
        if dmin<dist_min % only try to merge points that are not far apart
            w1=Weights(iclose); w2=Weights(iclose+1);
            R1=Radii(iclose); R2=Radii(iclose+1);
                % optimised barycenter
                A_Phik_GS=@(alpha,t,r) A_Phik(t,r)-w1*A_delta_R(t,R1,r)-w2*A_delta_R(t,R2,r) + (w1+w2)*A_delta_R(t,(1-alpha)*R1+alpha*R2,r);
                EDk_GS =@(alpha) s*integral2(@(t,r) t.^(s-1).*(1-A_Phik_GS(alpha,t,r)).^n .*psi(r,d),0,Inf,0,UB); % E{distortion}
                [aN,bN]=GoldenS(EDk_GS,0,1,1e-6); 
                alpha_star=(aN+bN)/2;
                A_Phik_try=@(t,r) A_Phik_GS(alpha_star,t,r);
                EDk_try=EDk_GS(alpha_star);
                Rnew=(1-alpha_star)*R1+alpha_star*R2;

            if EDk_try<EDk
                % merge
                A_Phik=@(t,r) A_Phik_try(t,r);
                EDk = EDk_try;
                Radii([iclose iclose+1])=[];
                Indices([iclose iclose+1])=[];
                Weights([iclose iclose+1])=[];
                RHO=[RHO Rnew]; % we have a new "grid point"
                LRHO=LRHO+1;
                Radii=[Radii Rnew];
                Indices=[Indices LRHO];
                Weights=[Weights w1+w2];
            end
        end
    end

    % update E distortion and sensitivity
    improvement=EDkold-EDk;
    EDkold=EDk;
    % sensitivity
    S_Phik=@(rho) n*s*integral2(@(t,r) t.^(s-1).*(1-A_Phik(t,r)).^(n-1) .*A_delta_R(t,rho,r) .*psi(r,d),0,Inf,0,UB);
    % check optimality
    Sens=NaN(1,LRHO);
    for i=1:LRHO
        rho=RHO(i);
        Sens(i)=S_Phik(rho);
    end
        
        [Radii
        Weights]
        improvement
        EDk
        Sens(Indices)
    S_bar_Phik=Weights*Sens(Indices)';
    one_minusEff=(max(Sens)-S_bar_Phik)/EDk
end    
ED=EDk;

% Remove negligible weights (since k_remove may be different from 1, they are not necessarily discarded at every iteration)
[~,i_negligible]=find(Weights<1e-12);
Radii(i_negligible)=[];
Indices(i_negligible)=[];
Weights(i_negligible)=[];
Weights=Weights/sum(Weights);

Sensitivities_on_RHO=Sens(1:LRHO_0);
Sensitivities_on_support=Sens(Indices);
end