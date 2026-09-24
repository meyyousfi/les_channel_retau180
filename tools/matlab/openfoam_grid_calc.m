function [X,Y,Z,x,y,z]=openfoam_grid_calc(sim,H,W,L,utau,nu,d1_plus,nx,ny,nz,span_mode,ar_top_surface)

% calculate input for blockMeshDict for flume or periodic channel

% input data
switch sim
    case 'channel'
    case 'flume'
%         H=0.1014;
%         W=1;
%         L=20*H;
%         utau=0.02119;
%         d1_plus=50;
%         span_mode='periodic';
%         % span_mode='wall';
%         ar_top_surface=3; % aspect ratio in the wall-normal direction - still square in span and streamwise
end


switch sim
    case 'channel'

        x=linspace(0,L,nx);
        disp(['dx uniform= ',num2str(x(2)-x(1))])
        disp(['dx+ uniform= ',num2str((x(2)-x(1))*utau/nu)])
        
        % solve the series (y_j+1-y_j)=alpha(y_j-y_j-1)-> y_j+1-y_j=a^(i-1)y_2
        % y is the coordinate of the cell faces (y_1=0 y_nj+1=H)
        % H/y_2=(alpha^nj-1)/(alpha-1) -> solve for alpha
        dy1=2*d1_plus*nu/utau; % cell face value - height of the first cell
        gr = fzero(@(a) (a^(ny-1)-1)/(a-1)-H/dy1,1.2); % growth ratio
        disp(['effective growth ratio= ',num2str(gr)])
        y=zeros(1,ny);
        y(2)=dy1;
        for j=2:ny-1
            y(j+1)=y(j)+gr*(y(j)-y(j-1));
        end
        % max factor for simplegrading
        fy=diff(y);
        disp(['dyplus at wall= ',num2str(fy(1)*utau/nu)])
        disp(['dyplus at top= ',num2str(fy(end)*utau/nu)])
        fy=fy(end)/fy(1);
        disp(['Simplegrading Factor Y= ',num2str(fy)])
        disp('*******************************************************************')
        
        z=linspace(0,W,nz);
        disp(['dz uniform= ',num2str(z(2)-z(1))])
        disp(['dz+ uniform= ',num2str((z(2)-z(1))*utau/nu)])
        
        [X,Y,Z]=ndgrid(x,y,z);
    case 'flume'        
        %% Wall-normal direction - condition y+<1 and growth ratio 1.04
        gr=1.04;
        dz1=d1_plus*nu/utau;
        
        % H*(gr-gr^2)-dz1*(gr-gr^Nz)=0 solution of series H=dz1*sum_(i=1)^(Nz-1)gr^(i-1)
        
        nz=ceil(log(-(H*(gr-gr^2)-gr*dz1)/dz1)/log(gr));
        % take nearest power of 2
        nz=2^nextpow2(nz); % # of cell centers
        disp(['Nz input= ',num2str(nz)])
        nz=nz+1; % # of cell faces
        
        % updated gr (solve sereis backward)
        gr = fzero(@(gr) H*(gr-gr^2)-dz1*(gr-gr^nz),gr);
        disp(['effective growth ratio= ',num2str(gr)])
        
        % check
        z=zeros(1,nz);
        z(1)=0;
        z(2)=dz1;
        for i=3:nz
            dz=z(i-1)-z(i-2);
            z(i)=z(i-1)+gr*dz;
        end
        % disp(['z= ',num2str(z)])
        % disp(['len(z)= ',num2str(numel(z))])
        
        % max factor for simplegrading
        fz=diff(z);
        disp(['dzplus at wall= ',num2str(fz(1)*utau/nu)])
        disp(['dzplus at top= ',num2str(fz(end)*utau/nu)])
        fz=fz(end)/fz(1);
        disp(['Simplegrading Factor Z= ',num2str(fz),' for the only block 1'])
        disp('*******************************************************************')
        
        %% streamwise direction - condition: square cells at the top boundary
        dy=ar_top_surface*(z(nz)-z(nz-1));
        ny=L/2/dy;
        % take nearest power of 2
        ny=2^nextpow2(ny); % # of cell centers
        disp(['ny= ',num2str(ny),' input per block'])
        ny=ny+1; % # of cell faces
        
        % check
        y=linspace(0,L,ny*2-1);
        % disp(['y= ',num2str(y)])
        % disp(['len(y)= ',num2str(numel(y))])
        % max factor for simplegrading
        fy=diff(y);
        disp(['dyplus at inlet= ',num2str(fy(1)*utau/nu)])
        disp(['dyplus at outlet= ',num2str(fy(end)*utau/nu)])
        fy=fy(end)/fy(1);
        disp(['Simplegrading Factor Y= ',num2str(fy),' for both blocks'])
        disp('*******************************************************************')
        
        %% spanwise direction
        
        switch span_mode
            case 'periodic' % - condition: square cells at the top boundary
                dx=ar_top_surface*(z(nz)-z(nz-1));
                nx=W/2/dx;
                % take nearest power of 2
                nx=2^nextpow2(nx); % # of cell centers
                disp(['nx= ',num2str(nx),' input per block'])
                nx=nx+1; % # of cell faces
                
                % check
                x=linspace(0,W,nx*2-1);
                fx=diff(x);
                disp(['dxplus at side= ',num2str(fx(1)*utau/nu)])
                disp(['dxplus at channel center= ',num2str(fx(nx)*utau/nu)])
                fx=fx(end)/fx(1);
                disp(['Simplegrading Factor X= ',num2str(fx),' for block 1 and 2'])
            case 'wall' % - condition: condition y+<1 and growth ratio 1.04
                gr=1.04;
                dx1=dz1;
                
                % W/2*(gr-gr^2)-dx1*(gr-gr^Nx)=0 solution of series W/2=dx1*sum_(i=1)^(Nx-1)gr^(i-1)
                
                nx=ceil(log(-(W/2*(gr-gr^2)-gr*dx1)/dx1)/log(gr));
                % take nearest power of 2
                nx=2^nextpow2(nx); % # of cell centers
                disp(['Nx= ',num2str(nx),' input per block'])
                nx=nx+1; % # of cell faces
                
                % updated gr (solve sereis backward)
                gr = fzero(@(gr) W/2*(gr-gr^2)-dx1*(gr-gr^nx),gr);
                disp(['effective growth ratio= ',num2str(gr)])
                % check
                x=zeros(1,nx*2-1);
                x(2)=dx1;
                for i=3:nx
                    dx=x(i-1)-x(i-2);
                    x(i)=x(i-1)+gr*dx;
                end
                x(nx*2-1)=W;
                x(nx*2-2)=W-dx1;
                for i=nx*2-3:-1:nx+1
                    dx=x(i+2)-x(i+1);
                    x(i)=x(i+1)-dx*gr;
                end
                fx=diff(x);
                disp(['dxplus at side= ',num2str(fx(1)*utau/nu)])
                disp(['dxplus at channel center= ',num2str(fx(nx)*utau/nu)])
                fx=fx(nx)/fx(1);
                disp(['Simplegrading Factor X= ',num2str(fx),' for block 1'])
                disp(['Simplegrading Factor X= ',num2str(1/fx),' for block 2'])
        end
        % max factor for simplegrading
        % disp(['x= ',num2str(x)])
        % disp(['len(x)= ',num2str(numel(x))])
        
        
        disp(['Ntot (total # of cells)= ',num2str((numel(x)-1)*(numel(y)-1)*(numel(z)-1)/1.e6),' mln points'])
end

% plot_grid(x,y,z);
end

function plot_grid(x,y,z)
close all

figure
[X,Y]=meshgrid(x,y);
plot(X,Y,'k',X',Y','k','LineWidth',1)
xlabel 'x'
ylabel 'y'
axis equal

figure
[X,Z]=meshgrid(x,z);
plot(X,Z,'k',X',Z','k','LineWidth',1)
xlabel 'x'
ylabel 'z'
axis equal

figure
[Y,Z]=meshgrid(y,z);
plot(Y,Z,'k',Y',Z','k','LineWidth',1)
xlabel 'y'
ylabel 'z'
axis equal
end

