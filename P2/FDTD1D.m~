%FDTD1D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Pre-Program Work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize MATLAB
close all; clc;
clear all; 

%Constants
c0 = 299792458; %m/s
e0 = 8.854187817*10^-12; %F/m
u0 = 1.256637061*10^-6; %H/m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Simulated Environment Settings
STEPS = 1000;
f_max = 1e9;  % 1Ghz
Nz = 180;
dz = 1.4286e-8;

%Compute Time Steps
dt = dz/(2*c0); %secs

%Grid Axis
za=[0:Nz-1]*dz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Material Vectors
ER = ones([1 Nz]);
UR = ones([1 Nz]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nzc = round (Nz/2);  %Position of Sources
tau = 1.2e-15        %normally 0.5/f_max *** Not working***
ta = [0:STEPS-1]*dt;     % Time Axis;
t0 = 6*tau;              % Delay
s = dz/(2*c0) + dt/2;    % Delay between E and H
Esrc = exp(-((ta-t0)/tau).^2); % E Source
A = -sqrt(ER(nzc)/UR(nzc));    % H Amplitude
Hsrc = A*exp(-((ta-t0+s)/tau).^2); % H Source

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FDTD Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Update Coefficients
mER = (c0*dt/dz)./ER;
mHR = (c0*dt/dz)./UR;

% Initialize Feilds
Ey = zeros([1 Nz]);
Hx = zeros([1 Nz]);


%PAB Parameters
h1 = 0; h2 = 0; h3 = 0;
e1 = 0; e2 = 0; e3 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:STEPS
 
  % Calculate H
  for nz = 1:Nz-1
    Hx(nz) = Hx(nz) + mHR(nz)*(Ey(nz+1)-Ey(nz));
  end
  
  Hx(Nz) = Hx(Nz) + mHR(Nz)*(e3 - Ey(Nz));

  h3 = h2; h2 = h1; h1 = Hx(1); % Boundary Params;
  
  % Calculate E  
  Ey(1) = Ey(1) + mER(1)*(Hx(1) - h3);
  for nz = 2:Nz
    Ey(nz) = Ey(nz) + mER(nz)*(Hx(nz)-Hx(nz-1)); 
  end
  
  %Inject Source
  Ey(nzc) = Ey(nzc) + Esrc(t);

  e3=e2; e2=e1; e1=Ey(Nz); % Boundary Params;
 
  h = plot(za, Ey, '-b'); hold on;
  plot(za, Hx, '-r'); hold off;
  axis([za(1) za(Nz) -1.1 1.1]);
  xlabel('z');
  title(['Field at Step ' num2str(t) ' of ' num2str(STEPS)]);    
  drawnow();
    
  if(mod(t,50) == 0)
    saveas(h, ['images/' num2str(t) '.jpg'], 'jpg');
  end
 end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
SetFigure(fig, 'HW#3-P2', [500 274 965 826]);

%Plot Magnetic Field
subplot(211)
h = plot(Hx, '-r', 'LineWidth', 2);
title('Magnetic Field');
h = get(h, 'Parent');
set(h, 'Fontsize', 14);
xlabel('z');
ylabel('Hx', 'Rotation', 0);
set(gca,'YTickLabel',{'1','0.5','0', '-0.5', '-1'})

%Plot Electric Field
subplot(212)
h = plot(Ey, '-b', 'LineWidth', 2);
title('Electric Field');
h = get(h, 'Parent');
set(h, 'Fontsize', 14);
xlabel('z');
ylabel('Ey', 'Rotation', 0);
set(gca,'YTickLabel',{'1','0.5','0', '-0.5', '-1'})




