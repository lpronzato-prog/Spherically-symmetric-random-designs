% Example of construction of an optimal design (mixture of uniform distributions on concentric spheres)

%--------------------------------------------------------------------------
clear variables
rng('default') % -> always the same random numbers
%rng('shuffle');
clear all
%--------------------------------------------------------------------------
% Colors
Goldenrod=[1 0.9 0.15];
red=[0.9,0.0,0.0];
yellow="#FAB84D";
blue="#3DA3FF";
green="#30E172";
purple="#8D4AC8";
grey1="#F2F2F2";
grey2="#C8C8C8";
grey1="#848484";
grey4="#585858";
%--------------------------------------------------------------------------
% Choose a distribution

%distribution='ball';    % mu uniform in B_d(0,1)
distribution='normal';  % mu normal N(0,Id/d)

sigma=@(t,R,r) min(max((t.^2-(-r+R).^2)./(4*r*R),0),1); % first parameter for beta inverse (truncated to [0,1])
switch lower(distribution)
    case 'ball'
        psi=@(r,d) d*r.^(d-1);                  % U uniform in B_d(0,1)
        UB=1;                                   % upper bound for integration
        UB_GS=1;                                % upper bound for GS algorithm (to find the optimal radius of one sphere)
        RHO_0=0:0.01:1;                           % set of radii considered 
    case 'normal'
        psi=@(r,d) 1/(2^(d/2-1))*d^(d/2)/gamma(d/2)*r.^(d-1).*exp(-d*r.^2/2); % U normal N(0,Id/d)
        UB=Inf;                                 % upper bound for integration
        UB_GS=2;                                % upper bound for GS algorithm (to find the optimal radius of one sphere)
        RHO_0=0:0.03:2;                         % set of radii considered 
end
LRHO_0=length(RHO_0);
one_minuP_to_n=@(t,R,r,d,n) ( 1-betainc(sigma(t,R,r),(d-1)/2,(d-1)/2) ).^n;

% P=Pa: X uniform on S_{d-1}(a);

ED=@(R,d,n,s) s*integral2(@(t,r) t.^(s-1).*one_minuP_to_n(t,R,r,d,n).*psi(r,d),0,Inf,0,UB); % E{distortion}

%__________________________________________________________________________
% illustrate the optimality of a_star by plotting the sensitivity function,
% compare with the sensitivity for an arbitrary value a_dum
a_dum=0.5; % an arbitrary sphere radius, not optimal

switch lower(distribution)
    case 'ball' % mu uniform in B_d(0,1)
        n=100
        %n=509  % 509 is the design size for which the uniform distribution
                % on a single sphere stops being optimal for s=2 and d=8
        s=2
        d=8   
        ED_a=@(R) ED(R,d,n,s);
        [aN,bN]=GoldenS(ED_a,1e-2,UB_GS,1e-3);        
    case 'normal' % mu normal N(0,Id/d)
        n=100
        %n=176  % 176 is the design size for which the uniform distribution
                % on a single sphere stops being optimal for s=5 and d=10
        s=5
        d=10
        ED_a=@(R) ED(R,d,n,s);
        [aN,bN]=GoldenS(ED_a,1e-2,UB_GS,1e-3);     
end
a_star=(aN+bN)/2 
ED_star=ED_a(a_star)

% Check optimality of this a_star
A_delta_R=@(t,R,r) betainc(sigma(t,R,r),(d-1)/2,(d-1)/2);

ED_bis=@(R) s*integral2(@(t,r) t.^(s-1).*(1-A_delta_R(t,R,r)).^n .*psi(r,d),0,Inf,0,UB); % E{distortion}
ED_bis(a_star) % coincides with ED

% compute sensitivity
S_delta_a=@(a,rho) n*s*integral2(@(t,r) t.^(s-1).*(1-A_delta_R(t,a,r)).^(n-1) .*A_delta_R(t,rho,r) .*psi(r,d),0,Inf,0,UB);

RHO=RHO_0;
LRHO=length(RHO);

Sens=NaN(1,LRHO);
Sens_star=NaN(1,LRHO);

for i=1:LRHO
    rho=RHO(i);
    Sens(i)=S_delta_a(a_dum,rho);
    Sens_star(i)=S_delta_a(a_star,rho);
end
Sens_bar_star=S_delta_a(a_star,a_star)
Sens_bar_dum=S_delta_a(a_dum,a_dum)

sensmin=min(min(Sens),min(Sens_star));
figure(1)
n0=10; step=10;
pp1=plot(RHO,Sens_star,'-',...
    'LineWidth',2,'MarkerSize',7,'MarkerIndices',1:step:LRHO-n0);
hold on
set(pp1,'Color',blue,'LineStyle','-','Marker','o','MarkerFaceColor',blue)
pp11=plot([0 max(RHO)],[Sens_bar_star Sens_bar_star],'k-','LineWidth',2);
set(pp11,'Color',blue,'LineStyle',':','Marker','o','MarkerFaceColor',blue)
pp12=plot([a_star a_star],[sensmin Sens_bar_star],'-','LineWidth',2);
set(pp12,'Color',blue,'LineStyle','--','Marker','o','MarkerFaceColor',blue)


pp2=plot(RHO,Sens,'-',...
    'LineWidth',2,'MarkerSize',7,'MarkerIndices',1:step:LRHO-n0);
set(pp2,'Color',yellow,'LineStyle','-','Marker','diamond','MarkerFaceColor',yellow)
pp21=plot([0 max(RHO)],[Sens_bar_dum Sens_bar_dum],'k-','LineWidth',2);
set(pp21,'Color',yellow,'LineStyle',':','Marker','diamond','MarkerFaceColor',yellow)
pp22=plot([a_dum a_dum],[sensmin Sens_bar_dum],'k--','LineWidth',2);
set(pp22,'Color',yellow,'LineStyle','--','Marker','diamond','MarkerFaceColor',yellow)

xlabel('$\rho$','fontsize',20,'Interpreter','latex')
set(gca,'FontSize',20)
axis tight
grid on
hold off

%__________________________________________________________________________
%% construct an optimal mixture
switch lower(distribution)
    case 'ball' % mu uniform in B_d(0,1)
        n=10000, epsi=1e-5;        
        s=2
        d=8
        initial='no';
    case 'normal' % mu normal N(0,Id/d)
        n=1000, epsi=1e-5;         
        s=5
        d=10
        initial='1pt-opt';
end

RHO=RHO_0;
LRHO=length(RHO);

[Radii, Weights,EDk,Sensitivities_on_support,Sensitivities_on_RHO] = Optimal_Sphere_Mixture(distribution,d,n,s,RHO,epsi,1,1,5,1e-3,1,0.2,initial);

%--------------------------------------------------------------------------
% Plot the cdf of Zador's asymptotically optimal distribution 
% and of the optimal, non-asymptotic, one

[rho_sort,isort]=sort(Radii,'ascend');
w_sort=Weights(isort);
FF=cumsum(w_sort);
[rho_sort
 w_sort]

x_emp=sort([Radii Radii]);
F_emp=sort([0 FF(1:end-1) FF]);

switch lower(distribution)
    case 'ball' % mu uniform in B_d(0,1)
        figure(2)
        pp3=plot(RHO_0,RHO_0.^d,'-',...
            'LineWidth',2,'MarkerSize',7,'MarkerIndices',1:step:LRHO_0-n0);
        set(pp3,'Color',yellow,'LineStyle','--','Marker','diamond','MarkerFaceColor',yellow)
        hold on
        pp1=plot(x_emp,F_emp,'-','LineWidth',2); 
        set(pp1,'Color',blue,'LineStyle','-')
        pp2=plot(rho_sort,FF,'.');
        set(pp2,'Color',blue,'Marker','o','MarkerSize',10,'MarkerFaceColor',blue)
        xlabel('$\rho$','fontsize',20,'Interpreter','latex')
        set(gca,'FontSize',20)
        axis tight
        grid on
        hold off    
    case 'normal' % mu normal N(0,Id/d)
        figure(2)
        sig=sqrt(1+s/d);
        Fasympt=@(R) integral(@(r) 1/sig*psi(r/sig,d),0,R);
        FFasympt=NaN(1,length(RHO_0));
        for i=1:length(RHO_0)
            rho=RHO_0(i);
            FFasympt(i)=Fasympt(rho);
        end
        pp3=plot(RHO_0,FFasympt,'-',...
            'LineWidth',2,'MarkerSize',7,'MarkerIndices',1:step:LRHO_0-n0);
        set(pp3,'Color',yellow,'LineStyle','--','Marker','diamond','MarkerFaceColor',yellow)
        hold on
        pp1=plot(x_emp,F_emp,'-','LineWidth',2); 
        set(pp1,'Color',blue,'LineStyle','-')
        pp2=plot(rho_sort,FF,'.');
        set(pp2,'Color',blue,'Marker','o','MarkerSize',10,'MarkerFaceColor',blue)
        xlabel('$\rho$','fontsize',20,'Interpreter','latex')
        set(gca,'FontSize',20)
        axis tight
        grid on
        hold off
end

%--------------------------------------------------------------------------
% Plot the sensitivity function and check optimality
figure(3)
n0=10; step=10;
pp1=plot(RHO_0,Sensitivities_on_RHO,'-',...
    'LineWidth',2);
hold on
set(pp1,'Color',blue,'LineStyle','-')

for i=1:length(Radii)
    plot([Radii(i) Radii(i)],[min(Sensitivities_on_RHO) Sensitivities_on_support(i)],'k--','LineWidth',2)
end

xlabel('$\rho$','fontsize',20,'Interpreter','latex')
set(gca,'FontSize',20)
%axis tight
grid on
hold off


% Compute the efficiency (in terms of expected distortion) of Zador's asymptotically optimal distribution  
switch lower(distribution)
    case 'ball'
        psiball=@(r,d) d*r.^(d-1);                  % U uniform in B_d(0,1)
        
        % integral_P_R_ball(tt,R,rr,d) computes 
        %   integral_P_R=@(t,R,r,d) integral(@(R) betainc(sigma(t,R,r),(d-1)/2,(d-1)/2).*phiball_design(R,d,sig),0,R);
        %   for all t=T(i) and r=rr(i), i=1,...,length(T)
        % E distortion
        R=1;
        ED_asympt_ball= s*integral2(@(t,r) t.^(s-1).*(1-integral_P_R_ball(t,R,r,d)).^n.*psiball(r,d),0,Inf,0,1) ; 
        Eff_asympt_ball = EDk/ED_asympt_ball 

    case 'normal'
        psinormal=@(r,d) 1/(2^(d/2-1))*d^(d/2)/gamma(d/2)*r.^(d-1).*exp(-d*r.^2/2);
        
        % integral_P_R_normal(tt,sig,rr,d) computes 
        %   integral_P_R=@(t,sig,r,d) integral(@(R) betainc(sigma(t,R,r),(d-1)/2,(d-1)/2).*phinormal_design(R,d,sig),0,Inf);
        %   for all t=T(i) and r=rr(i), i=1,...,length(T)
        % E distortion
        sig=sqrt(1+s/d); % From Zador
        ED_asympt_normal= s*integral2(@(t,r) t.^(s-1).*(1-integral_P_R_normal(t,sig,r,d)).^n.*psinormal(r,d),0,Inf,0,Inf) ;
        Eff_asympt_normal = EDk/ED_asympt_normal
end