Defs
PlotSpec

h=1;
W=pi;
L=3*pi;
ut=1;
ret=180;
nu=ut*h/ret;
d1_plus=0.15; % cell center value
nx=96;
ny=96;
nz=96;
[X,Y,Z]=openfoam_grid_calc('channel',h,W,L,ut,nu,d1_plus,nx,ny,nz,[],[]);

return

kappa=0.41;
b=5;
dp=11; % y+ location where log law and viscous law meet

%% check bulk
[ub]=calc_ub_smooth(kappa,b,ut,nu,h,dp);
reb=(ret/.09)^(1/0.88);
ub2=reb*nu/2/h;
disp(['Difference between ub from integral and from Pope law= ',num2str(abs(ub-ub2)/ub2*100),'%'])
disp(['ub= ',num2str(ub)])

%% build synthetic velocity profile
%% create mean profile
% Define transition range
y=Y(1,:,1);
Y_lower = 7*nu/ut;
Y_upper = 15*nu/ut;

% Two profile functions
f = tanh_smooth(y, Y_lower, Y_upper);
vel1 = y .* ut^2 / nu;
vel2 = zeros(size(y));
vel2(2:end) = ut / kappa * log(y(2:end) * ut / nu) + b * ut;

% Initialise output
vel_prof_th = zeros(size(y));
vel_prof_th = f.*vel1 + (1-f).*vel2;

%% add synthetic fluctuations (normal random distribution)
% value = min + (max - min) * rand;
% -1+2*rand generates random values between -1 and 1
% f = asymmetric_peak(Y(1,:,1), dp*nu/ut, max(Y(:))) generates a peak at dp*nu/ut
f = asymmetric_peak(Y(1,:,1), dp*nu/ut, max(Y(:)),2);
f=repmat(f,[nx,1,nz]);
u=zeros(size(X));
u=repmat(vel_prof_th,[nx,1,nz]);
min_val=-0.3*ub;
max_val= 0.3*ub;
ind=Y<dp*nu/ut;
u(ind)=u(ind) + (min_val + (max_val - min_val) * rand(size(X(ind)))).*f(ind);
min_val=-0.21*ub;
max_val= 0.21*ub;
ind=Y>=dp*nu/ut & Y<2*dp*nu/ut;
u(ind)=u(ind) + (min_val + (max_val - min_val) * rand(size(X(ind)))).*f(ind);
min_val=-0.12*ub;
max_val= 0.12*ub;
ind=Y>=2*dp*nu/ut;
u(ind)=max(u(ind)+ (min_val + (max_val - min_val) * rand(size(X(ind)))).*f(ind),0);
u(:,1,:)=0;
u(1,:,:)=u(end,:,:);
u(:,:,1)=u(:,:,end);
min_val=-0.1*ub;
max_val= 0.1*ub;
v=zeros(size(X)) + (min_val + (max_val - min_val) * rand(size(X))).*f;
v(:,1,:)=0;
v(1,:,:)=v(end,:,:);
v(:,:,1)=v(:,:,end);
min_val=-0.1*ub;
max_val= 0.1*ub;
w=zeros(size(X)) + (min_val + (max_val - min_val) * rand(size(X))).*f;
w(:,1,:)=0;
w(1,:,:)=w(end,:,:);
w(:,:,1)=w(:,:,end);

fprintf('Original v field average: %.4f\n', mean(v(:)));
fprintf('Original w field average: %.4f\n', mean(w(:)));

% Parameters
smooth_sigma = 0.2;      % Gaussian smoothing parameter (larger = more smoothing)

% Step 1: Smoothen the field to reduce fluctuations
fprintf('Smoothing the 3D field...\n');
u = smooth3(u, 'gaussian', [3 3 3], smooth_sigma);
v = smooth3(v, 'gaussian', [3 3 3], smooth_sigma);
w = smooth3(w, 'gaussian', [3 3 3], smooth_sigma);
fprintf('After smoothing v average: %.4f\n', mean(v(:)));
fprintf('After smoothing w average: %.4f\n', mean(w(:)));

% Alternative smoothing methods:
% field_smoothed = smooth3(field_3d, 'box', 5);  % Box filter
% field_smoothed = imgaussfilt3(field_3d, smooth_sigma);  % Requires Image Processing Toolbox

% Step 2: Normalize to achieve target average
target_average = 0.0;  % Desired average value for the field
fprintf('Normalizing the field...\n');
v = v - mean(v(:)) + target_average;
w = w - mean(w(:)) + target_average;

% Verify the result
fprintf('Final normalized v average: %.4f\n', mean(v(:)));
fprintf('Final normalized w average: %.4f\n', mean(w(:)));
fprintf('Target v and w average: %.4f\n', target_average);

% check rms of velocities
temp=u-repmat(vel_prof_th,[nx,1,nz]);
rms_u=sqrt(sum(temp.^2,'all')/(numel(u)-1));
fprintf('Fluctuation rms of u field (/ub): %.4f\n', rms_u/ub);
temp=v-target_average;
rms_v=sqrt(sum(temp.^2,'all')/(numel(v)-1));
fprintf('Fluctuation rms of v field (/ub): %.4f\n', rms_v/ub);
temp=w-target_average;
rms_w=sqrt(sum(temp.^2,'all')/(numel(w)-1));
fprintf('Fluctuation rms of w field (/ub): %.4f\n', rms_w/ub);

% plot fields
fig=2;
figure(fig)
subplot(3,1,1)
plot_2d(fig,X(:,:,1),Y(:,:,1),squeeze(u(:,:,1)),'x','y','U',false);
subplot(3,1,2)
plot_2d(fig,X(:,:,1),Y(:,:,1),squeeze(v(:,:,1)),'x','y','V',false);
subplot(3,1,3)
plot_2d(fig,X(:,:,1),Y(:,:,1),squeeze(w(:,:,1)),'x','y','W',false);
[cb]=set_contour_levels_wholefigure(gcf,'vel','best','','average','vik');

fig=3;
figure(fig)
n_plus=squeeze(Y(1,:,1))*ut/nu;
u_plus=squeeze(u(1,:,1))/ut;
% u_plus=squeeze(mean(u,[1 3]))/ut;
plot_1d(fig,n_plus,u_plus,'semilogx','','','-','o','b','Generated initial profile',2,1,'best','on',6);
plot_1d(fig,n_plus,squeeze(vel_prof_th(1,:,1)),'semilogx','','','-','','k','Mean theoreticalprofile',2,1,'best','on',[]);

n_plus_theory = logspace(-1, 3, 100);

% Viscous sublayer: u+ = n+
viscous_layer = n_plus_theory;

% Log layer: u+ = 1/kappa * ln(n+) + C
kappa = 0.41;  % von Karman constant
C = 5.0;       % Constant
log_layer = (1/kappa) * log(n_plus_theory) + C;

% Plot theoretical lines
plot_1d(fig,n_plus_theory,viscous_layer,'semilogx','','','--','','r',...
    'u^+ = n^+ (Viscous sublayer)',1.5,1,'','on',[]);
plot_1d(fig,n_plus_theory,log_layer,'semilogx','','','--','','g',...
    sprintf('u^+ = %.2f ln(n^+) + %.1f (Log layer)', 1/kappa, C),1.5,1,'best','on',[]);
grid on;
xlabel('n^+ (Wall distance in wall units)', 'FontSize', 12, 'Interpreter', 'tex');
ylabel('u^+ (Velocity in wall units)', 'FontSize', 12, 'Interpreter', 'tex');
xlim([0.1, max(n_plus)*1.1]);
ylim([0, max(u_plus)*1.1]);

%% velocity variance
u_prime=u-repmat(mean(u,[1,3]),[nx,1,nz]);
v_prime=v-repmat(mean(v,[1,3]),[nx,1,nz]);
w_prime=w-repmat(mean(w,[1,3]),[nx,1,nz]);
uu_prime=u_prime.^2;
vv_prime=v_prime.^2;
ww_prime=w_prime.^2;

fig=4;
figure(fig);
phi=squeeze(mean(uu_prime,[1 3]))/ut^2;
plot_1d(fig,n_plus,phi,'semilogx','','','-','',LineColor(1,:),...
    '$u\prime u\prime/u_\tau^2$',1.5,1,'','on',[]);
phi=squeeze(mean(vv_prime,[1 3]))/ut^2;
plot_1d(fig,n_plus,phi,'semilogx','','','-','',LineColor(2,:),...
    '$v\prime v\prime/u_\tau^2$',1.5,1,'','on',[]);
phi=squeeze(mean(ww_prime,[1 3]))/ut^2;
plot_1d(fig,n_plus,phi,'semilogx','$n^+$','','-','',LineColor(3,:),...
    '$w\prime w\prime/u_\tau^2$',1.5,1,"best",'on',[]);

%% write file for openFOAM input
% Define the filename
filename = 'initialVelocityProfile';

% Define the number of vectors
numVectors = (nz)*(ny)*(nz);

% Open the file for writing
fid = fopen(filename, 'w');

% Write the header
fprintf(fid, 'velocityProfile\n');
fprintf(fid, 'nonuniform List<vector>\n');
fprintf(fid, '%d\n', numVectors);
fprintf(fid, '(\n');

% Write each velocity vector
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            fprintf(fid, '(%f %f %f)\n', u(i, j, k), v(i, j, k), w(i, j, k));
        end
    end
end

% Write the closing bracket
fprintf(fid, ');\n');

% Close the file
fclose(fid);

disp(['Velocity File written to ', filename]);

%% pressure field
fprintf('Computing pressure field...\n');
%    F_components   Cell array of scalar field arrays, one per dimension:
%                     2-D : {Fx [Nx×Ny],  Fy [Nx×Ny]}
%                     3-D : {Fx [Nx×Ny×Nz], Fy [Nx×Ny×Nz], Fz [Nx×Ny×Nz]}
%    coords         Cell array of strictly-increasing node coordinate vectors:
%                     2-D : {x [Nx×1],  y [Ny×1]}
%                     3-D : {x [Nx×1],  y [Ny×1],  z [Nz×1]}
%    bc             Struct with per-face boundary condition strings.
%                   Required fields:
%                     2-D : x_low, x_high, y_low, y_high
%                     3-D : x_low, x_high, y_low, y_high, z_low, z_high
%    bc = struct('x_low','periodic','x_high','periodic', ...
%                'y_low','periodic','y_high','periodic', ...
%                'z_low','periodic','z_high','periodic');

rho=1;
fprintf('Setting fluid density: %.4f\n', rho);
delta_t=0.002;
fprintf('Setting time step: %.4f\n', delta_t);
dp_dx=-rho*ut^2/h;
fprintf('Calculating pressure gradient: %.4f\n', dp_dx);
x=X(:,1,1);
y=Y(1,:,1);
z=Z(1,1,:);
coords={x,y,z};
u_bc = struct('x_low','periodic','x_high','periodic', ...
              'y_low','wall','y_high','symmetry', ...
              'z_low','periodic','z_high','periodic');
divU = compute_divergence({u,v,w}, coords, u_bc);
p_bc = struct('x_low','periodic','x_high','periodic', ...
              'y_low','neumann','y_high','symmetry', ...
              'z_low','periodic','z_high','periodic');
f=rho/delta_t*divU;
p = solve_poisson(f, coords, p_bc);

p=p+dp_dx.*X-dp_dx*max(x(:));

fig=5;
figure(fig)
plot_2d(gcf,X(:,:,1),Y(:,:,1),squeeze(p(:,:,1)),'x','y','p',false);
[cb]=set_contour_levels_wholefigure(gcf,'$p$','best','','average','vik');
%% write file for openFOAM input
% Define the filename
filename = 'initialPressureProfile';

% Define the number of scalars
numScalars = (nz)*(ny)*(nz);

% Open the file for writing
fid = fopen(filename, 'w');

% Write the header
fprintf(fid, 'pressureProfile\n');
fprintf(fid, 'nonuniform List<scalar>\n');
fprintf(fid, '%d\n', numScalars);
fprintf(fid, '(\n');

% Write each velocity vector
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            fprintf(fid, '%f \n', p(i, j, k));
        end
    end
end

% Write the closing bracket
fprintf(fid, ');\n');

% Close the file
fclose(fid);

disp(['Pressure File written to ', filename]);

%% pressure field
fprintf('Computing sgs viscosity field...\n');

viscos=2*ub*h/reb;

cs=0.01;
fprintf('Setting sgs viscosity constant: %.4f\n', cs);
x=X(:,1,1);
y=Y(1,:,1);
z=Z(1,1,:);
coords={x,y,z};
u_bc = struct('x_low','periodic','x_high','periodic', ...
              'y_low','wall','y_high','symmetry', ...
              'z_low','periodic','z_high','periodic');
gradU = compute_jacobian({u,v,w}, coords, u_bc);
gradU(1,:,:)=gradU(end,:,:);
gradU(:,:,1)=gradU(:,:,end);
nd = size(gradU, 2);   % number of spatial dimensions
S  = cell(nd, nd);
S_norm = zeros(size(u));

for i = 1:nd
    for j = 1:nd
        S{i,j} = 0.5 * (gradU{i,j} + gradU{j,i});
        S_norm = S_norm + S{i,j}.^2;
    end
end
S_norm = sqrt(2*S_norm);
DV = compute_cell_volumes(coords);
DV=DV.^(1/3);
ls  = min(kappa*Y,cs*DV);
nut = ls.^2.*S_norm;

fig=6;
figure(fig)
f=squeeze(mean(nut,3))/viscos;
plot_2d(gcf,X(:,:,1),Y(:,:,1),f,'x','y','$\nu_{sgs}/\nu$',false);
[cb]=set_contour_levels_wholefigure(gcf,'$\nu_{sgs}/\nu$','best','','average','vik');

fig=7;
figure(fig)
f=squeeze(mean(nut,[1,3]))/viscos;
plot_1d(gcf,y,f,'plot','y','$\nu_{sgs}/\nu$','-','','k','Generated initial profile',2,1,'best','on',6);
%% write file for openFOAM input
% Define the filename
filename = 'initialViscosityProfile';

% Define the number of scalars
numScalars = (nz)*(ny)*(nz);

% Open the file for writing
fid = fopen(filename, 'w');

% Write the header
fprintf(fid, 'nutProfile\n');
fprintf(fid, 'nonuniform List<scalar>\n');
fprintf(fid, '%d\n', numScalars);
fprintf(fid, '(\n');

% Write each velocity vector
for k = 1:nz
    for j = 1:ny
        for i = 1:nx
            fprintf(fid, '%f \n', nut(i, j, k));
        end
    end
end

% Write the closing bracket
fprintf(fid, ');\n');

% Close the file
fclose(fid);

disp(['Viscosity File written to ', filename]);

%%


function [ub]=calc_ub_smooth(kappa,B,ut,nu,h,d1p)
    d1=d1p*nu/ut;
    re_t_s=h*ut/nu;
    ub=nu/(h-d1)*(re_t_s/kappa*(log(re_t_s)-1)+B*re_t_s)-...
       nu/(h-d1)*(d1p/kappa*(log(d1p)-1)+B*d1p);
end

function f = tanh_smooth(Y, y_lower, y_upper)
    % tanh-based smoothing function: f=1 at y_lower, f=0 at y_upper
    % Sharper transitions for smaller width
    width = (y_upper - y_lower) / 12; % Adjust divisor for sharpness
    midpoint = (y_lower + y_upper) / 2;
    f = 0.5 * (1 - tanh((Y - midpoint) / width));
end
