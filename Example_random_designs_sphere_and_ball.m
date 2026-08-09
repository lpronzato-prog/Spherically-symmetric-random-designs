% Example of construction of spherically symmetric random (or quasi-random) designs 
% in high dimension with small Ls-mean quantisation error

d=30; 
n=10000;
s=2;
accuracy=1e-2;

% The reference measure is uniform on the sphere S_{d-1}(0,1)
    disp('mu uniform on S_{d-1}(0,1), design uniform on S_{d-1}(0,R)')
    % For n random points uniformly distributed on S_{d-1}(0,R), the optimal radius is:
        Ropt = Ropt_sphere_sphere(d,n,s,accuracy)
    % the expected distortion (distortion = (Ls-mean quantisation error)^s) is:
        Edist = Edistortion_sphere_sphere(Ropt,d,n,s)

% The reference measure is uniform in the ball B_d(0,1)
    disp('mu uniform in B_d(0,1), design uniform in B_d(0,R)')
    % For n random points uniformly distributed in B_d(0,R), the optimal radius is:
        Ropt = Ropt_ball('ball',d,n,s,accuracy)
    % the expected distortion (distortion = (Ls-mean quantisation error)^s) is:
        Edist = Edistortion_ball_ball(Ropt,d,n,s)
    %
    disp('mu uniform in B_d(0,1), design uniform on S_{d-1}(0,R)')
    % For n random points uniformly distributed on S_{d-1}(0,R), the optimal radius is:
        Ropt = Ropt_ball('sphere',d,n,s,accuracy)
    % the expected distortion (distortion = (Ls-mean quantisation error)^s) is:
        Edist = Edistortion_ball_sphere(Ropt,d,n,s)

% The reference measure is normal N(0,eye(d)/d).
    disp('mu normal N(0,eye(d)/d), design normal N(0,sig^2*eye(d)/d)')
    % For n random points independent normal N(0,sig^2*eye(d)/d), the optimal std is:
        sigopt = Ropt_normal('normal',d,n,s,accuracy)
    % the expected distortion (distortion = (Ls-mean quantisation error)^s) is:
        Edist = Edistortion_normal_normal(sigopt,d,n,s)
    %
    disp('mu normal N(0,eye(d)/d), design uniform on S_{d-1}(0,R)')
    % For n random points uniformly distributed on S_{d-1}(0,R), the optimal radius is:
        Ropt = Ropt_normal('sphere',d,n,s,accuracy)
    % the expected distortion (distortion = (Ls-mean quantisation error)^s) is:
        Edist = Edistortion_normal_sphere(Ropt,d,n,s)
        