
%------------------ 50 Ohm Microstrip transmission line with ABC & Luebber's Tapered Source ----------------------
%------------------ Feed point, ABC and substrate coincide at x-z plane -----------------
%------------------ To plot waves on transmitting & receiving ends -----------------

clc;
close all;
clear all;

%----- Substrate and line dimensions for 50 Ohm----------
sub_l = 40e-3;
sub_h = 1.57e-3;
sub_w = 40e-3;

strip_l = sub_l;
strip_w = 4.66e-3;

%----- Substrate material and surrounding air ----------
eps0 = 8.854e-12;
meu0 = 4*pi*1e-7;
epsr = 2.3;                   % RT/Duroid material
meur = 1;
epss = eps0*epsr;             % substrate
meus = meu0*meur;
epsa = eps0;                  % air
meua = meu0;
epse = (epsa + epss)/2;       % Average eps

%---- EM wave -----------
c = 3e8;
vmin = c/sqrt(epsr*meur);
vmax = c;
fmax = 15e9; 
lamda = vmin/fmax;

%------ FDTD cell length and time step---------
%------- Substrate thickness is in Y-direction
dx = 0.5e-3;       % 0.5e-3 and 80 cells
dy = dx/2;
dz = dx; 
dt = 1/(vmax*sqrt(1/dx^2 + 1/dy^2 + 1/dz^2));

xoff = 0;
yoff = 0;
zoff = 0; 

Nx = round(sub_w/dx) + 2*xoff;
Ny = 34 + 2*yoff;                        % upward | 30 cells above the substrate
Nz = round(sub_l/dz) + 2*zoff;

%----------- Start and end coordinates of the subsrate 
sub_x1 = 1 + xoff;
sub_x2 = 80 + xoff;
sub_y1 = 1 + yoff; 				   % thickness starts from index 1
sub_y2 = 4 + yoff;                 % 4 cells in height
sub_z1 = 1 + zoff;
sub_z2 = 80 + zoff;

%----------- Start and end coordinates of the strip
strip_x1 = round(0.5*(Nx-strip_w/dx)) + 1;          % +1 for odd no. of cells
Nw = round(strip_w/dx);
strip_x2 = strip_x1 + Nw - 1;
strip_y1 = sub_y1;
strip_y2 = sub_y2 + 1; 			% strip is on the top surface of the substrate
strip_z1 = sub_z1;              % strip is starting at ABC
strip_z2 = sub_z2; 
y_mid = strip_y1 + round((strip_y2-strip_y1)/2);

%-------- Feed point of the strip
feed_x = strip_x1 + round(0.5*strip_w/dx) - 1; 		% center of strip
feed_y = strip_y1;
feed_z1 = strip_z1;                             % feed starts at the strip edge
feed_z2 = strip_z2;

%------- Initialise E,H,multipliers and ABC arrays ------------
ex = single(zeros(Nx,Ny,Nz));
ey = single(zeros(Nx,Ny,Nz));
ez = single(zeros(Nx,Ny,Nz));

hx = single(zeros(Nx,Ny,Nz));
hy = single(zeros(Nx,Ny,Nz));
hz = single(zeros(Nx,Ny,Nz));

const_ex = single(zeros(Nx,Ny,Nz));
const_ey = single(zeros(Nx,Ny,Nz));
const_ez = single(zeros(Nx,Ny,Nz));

const_hx = single(zeros(Nx,Ny,Nz));
const_hy = single(zeros(Nx,Ny,Nz));
const_hz = single(zeros(Nx,Ny,Nz));

%------- Default constants in update eqn
const_ex(1:Nx,1:Ny,1:Nz) = dt/(epsa*dx); 
const_ey(1:Nx,1:Ny,1:Nz) = dt/(epsa*dy);
const_ez(1:Nx,1:Ny,1:Nz) = dt/(epsa*dz);

%--------------------------------------
const_hx(1:Nx,1:Ny,1:Nz) = dt/(meua*dx);
const_hy(1:Nx,1:Ny,1:Nz) = dt/(meua*dy);
const_hz(1:Nx,1:Ny,1:Nz) = dt/(meua*dz);

%----------- Substrate---------
const_ex(sub_x1:sub_x2,sub_y1:sub_y2,sub_z1:sub_z2) = dt/(epss*dx); 
const_ey(sub_x1:sub_x2,sub_y1:sub_y2,sub_z1:sub_z2) = dt/(epss*dy);
const_ez(sub_x1:sub_x2,sub_y1:sub_y2,sub_z1:sub_z2) = dt/(epss*dz);

%----------- Top substrate-air interface has eps_avg
const_ex(sub_x1:sub_x2,strip_y2,sub_z1:sub_z2) = dt/(epse*dx);
const_ez(sub_x1:sub_x2,strip_y2,sub_z1:sub_z2) = dt/(epse*dz);


%----------- Metal boundaries have Etan=0 ---------
%------- Constants in E-update eqn
gx(1:Nx,1:Ny,1:Nz) = 1; gy=gx; gz=gy; 		% set default = 1 for non-metal areas
fx=gx; fy=gy;

%---- GND plane
gx(sub_x1:sub_x2,sub_y1,sub_z1:sub_z2) = 0;
gz(sub_x1:sub_x2,sub_y1,sub_z1:sub_z2) = 0;

%---- Stripline
gx(strip_x1:strip_x2,strip_y2,strip_z1:strip_z2) = 0;
gz(strip_x1:strip_x2,strip_y2,strip_z1:strip_z2) = 0;

%---- Tapered Source
for ind = sub_y1-1:sub_y2-2
    gx((feed_x-ind):(feed_x+ind),ind+2,feed_z1) = 0;
    gz((feed_x-ind):(feed_x+ind),ind+2,feed_z1) = 0;
end

%------ ABC constants
ex_abc = ex; ey_abc = ex_abc; ez_abc = ey_abc;
constz_abc = (vmax*dt-dz)/(vmax*dt+dz);
consty_abc = (vmax*dt-dy)/(vmax*dt+dy);
constx_abc = (vmax*dt-dx)/(vmax*dt+dx);

%------ Plot gx at the surface of the substrate to verify the model
f1 = figure('Name','gx Plot at the Surface');

figure(f1);
arr = gx(:,strip_y2,:);
surf(1:Nz,1:Nx,arr(1:Nx,1:Nz));

%------ Gaussian pulse ------------
Ts = 10*dt; 		% pulse width
t0 = 3*Ts;          % delay
N_steps = 1;		% default iteration steps
count = 0;

receive_Voltage = zeros(1e6,1);
transmit_Voltage = zeros(1e6,1);
receive_Current = zeros(1e6,1);
transmit_Current = zeros(1e6,1);

Zs = 50;         % Source Impedance in Ohms

%************** Iteration loop ****************
while(N_steps>0)  
    
    N_steps = input('\n\n Time steps(0 to quit): ');
    
    for n=1:N_steps         
        time = count*dt;
        count = count+1;
        clc;
        fprintf('\n\tTime step : %d',count);
        clc;
        fprintf('\n\tTime step : %d',count);       
       
        % Store values of voltage (integral of ey along y direction) for plottting outside the for loop
        transmit_Voltage(count) = (sum(ey(feed_x,strip_y1:strip_y2,feed_z1+5)))*dy;         % Sampled 5 cells inwards from feed points
        receive_Voltage(count) = (sum(ey(feed_x,strip_y1:strip_y2,feed_z2-5)))*dy;          % Note: feed_z1+5 and feed_z2-5

        %--------------- Ampere's Law at x = feed_x, y = mid -----------------
        transmit_Current(count)	= -((hx(feed_x,y_mid,feed_z1+5)-hx(feed_x,y_mid,feed_z1+6))*dx + (hz(feed_x,y_mid,feed_z1+6)-hz(feed_x-1,y_mid,feed_z1+6))*dz);
        receive_Current(count)	= -((hx(feed_x,y_mid,feed_z2-6)-hx(feed_x,y_mid,feed_z2-5))*dx + (hz(feed_x,y_mid,feed_z2-5)-hz(feed_x-1,y_mid,feed_z2-5))*dz);

        % ----- Sinusoidal signal of 5 GHz---------
        pulse = (exp(-((time-t0)^2)/(Ts^2)))*(sin(2*pi*25e9*time));
        ey(feed_x,feed_y,feed_z1) = (pulse/dy) + ((transmit_Current(count)*Zs)/dy); 
        
        %------------------------ compute H -------------------------
        j=1:Ny-1; k=1:Nz-1;
        hx(:,j,k) = hx(:,j,k)+fx(:,j,k).*(const_hz(:,j,k).*(ey(:,j,k+1)-ey(:,j,k))-const_hy(:,j,k).*(ez(:,j+1,k)-ez(:,j,k)));     
        i=1:Nx-1; k=1:Nz-1;
        hy(i,:,k) = hy(i,:,k)+fy(i,:,k).*(const_hx(i,:,k).*(ez(i+1,:,k)-ez(i,:,k))-const_hz(i,:,k).*(ex(i,:,k+1)-ex(i,:,k))); 
        i=1:Nx-1; j=1:Ny-1;
        hz(i,j,:) = hz(i,j,:)+const_hy(i,j,:).*(ex(i,j+1,:)-ex(i,j,:))-const_hx(i,j,:).*(ey(i+1,j,:)-ey(i,j,:));
        
        %------------------------- compute E ------------------------
        j=2:Ny-1; k=2:Nz-1;
        ex(:,j,k) = ex(:,j,k)+gx(:,j,k).*(const_ey(:,j,k).*(hz(:,j,k)-hz(:,j-1,k))-const_ez(:,j,k).*(hy(:,j,k)-hy(:,j,k-1))); 
        i=2:Nx-1; k=2:Nz-1;
        ey(i,:,k) = ey(i,:,k)+gy(i,:,k).*(const_ez(i,:,k).*(hx(i,:,k)-hx(i,:,k-1))-const_ex(i,:,k).*(hz(i,:,k)-hz(i-1,:,k)));  
        i=2:Nx-1; j=2:Ny-1;
        ez(i,j,:) = ez(i,j,:)+gz(i,j,:).*(const_ex(i,j,:).*(hy(i,j,:)-hy(i-1,j,:))-const_ey(i,j,:).*(hx(i,j,:)-hx(i,j-1,:)));  
            
        %----------------- ABC---------------------  
        %---------- XY plane 
        ex(:,:,1) = ex_abc(:,:,2) + constz_abc*(ex(:,:,2)-ex(:,:,1)); 
        ey(:,:,1) = ey_abc(:,:,2) + constz_abc*(ey(:,:,2)-ey(:,:,1));        
                        
        ex(:,:,Nz) = ex_abc(:,:,Nz-1) + constz_abc*(ex(:,:,Nz-1)-ex(:,:,Nz)); 
        ey(:,:,Nz) = ey_abc(:,:,Nz-1) + constz_abc*(ey(:,:,Nz-1)-ey(:,:,Nz));      
        
        %---------- XZ plane       
        %ex(:,1,:) = ex_abc(:,2,:) + consty_abc*(ex(:,2,:)-ex(:,1,:)); 
        %ez(:,1,:) = ez_abc(:,2,:) + consty_abc*(ez(:,2,:)-ez(:,1,:)); 
        
        ex(:,Ny,:) = ex_abc(:,Ny-1,:) + consty_abc*(ex(:,Ny-1,:)-ex(:,Ny,:)); 
        ez(:,Ny,:) = ez_abc(:,Ny-1,:) + consty_abc*(ez(:,Ny-1,:)-ez(:,Ny,:)); 
          
        %---------- YZ plane
        ey(1,:,:) = ey_abc(2,:,:)+constx_abc*(ey(2,:,:)-ey(1,:,:)); 
        ez(1,:,:) = ez_abc(2,:,:)+constx_abc*(ez(2,:,:)-ez(1,:,:));        
                      
        ey(Nx,:,:) = ey_abc(Nx-1,:,:)+constx_abc*(ey(Nx-1,:,:)-ey(Nx,:,:)); 
        ez(Nx,:,:) = ez_abc(Nx-1,:,:)+constx_abc*(ez(Nx-1,:,:)-ez(Nx,:,:));  
       
        ex_abc = ex;
        ey_abc = ey;  
        ez_abc = ez;  
       
    end  % for loop ends   
 
    f2 = figure('Name','Voltage');
    f3 = figure('Name','Current');
    
    %time_scale = (0:count-1)*dt;
    time_scale = 1:count;               % temporary
    
    % Plot Voltage vs time    
    figure(f2);
    subplot(1,2,1);
    plot(time_scale,transmit_Voltage(1:count),'-');
    title('Transmission End');
    grid on; grid minor;
    xlim([min(time_scale) max(time_scale)]);
    
    subplot(1,2,2);
    plot(time_scale,receive_Voltage(1:count),'-');
    title('Receiving End');
    grid on; grid minor;
    xlim([min(time_scale) max(time_scale)]);

    % Plot Current vs time
    figure(f3);
    subplot(1,2,1);
    plot(time_scale,transmit_Current(1:count),'-');
    title('Transmission End');
    grid on; grid minor;
    xlim([min(time_scale) max(time_scale)]);

    subplot(1,2,2);
    plot(time_scale,receive_Current(1:count),'-');
    title('Receiving End');
    grid on; grid minor;
    xlim([min(time_scale) max(time_scale)]);

end  % while loop ends       


% Compute Fourier Transform
V_fft_transmit = fft(transmit_Voltage(1:count));             
V_fft_receive = fft(receive_Voltage(1:count));   

I_fft_transmit = fft(transmit_Current(1:count));             
I_fft_receive = fft(receive_Current(1:count));  

Zin = abs(V_fft_transmit)./abs(I_fft_transmit);
Tin = (Zin-Zs)./(Zin+Zs);

f4 = figure('Name','FFT of Voltage');
f5 = figure('Name','FFT of Current');

freq = (0:count-1)./(count*dt);

figure(f4);
subplot(1,2,1);
plot(freq,abs(V_fft_transmit));
title('Transmission End');
grid on; grid minor;
xlim([min(freq) max(freq)/10]);
subplot(1,2,2);
plot(freq,abs(V_fft_receive));
title('Receiving End');
grid on; grid minor;
xlim([min(freq) max(freq)/10]);

figure(f5);
subplot(1,2,1);
plot(freq,abs(I_fft_transmit));
title('Transmission End');
grid on; grid minor;
xlim([min(freq) max(freq)/10]);
subplot(1,2,2);
plot(freq,abs(I_fft_receive));
title('Receiving End');
grid on; grid minor;
xlim([min(freq) max(freq)/10]);

f6 = figure('Name','Input Transmission Coeeficient (Tin) vs. Frequency');
figure(f6);
plot(freq,20*log10(abs(Tin)));
grid on; grid minor;
xlim([min(freq) max(freq)/4]);
title('Input Transmission Coeeficient (Tin) vs. Frequency');
xlabel('Frequency');
ylabel('20 log |Tin|');
