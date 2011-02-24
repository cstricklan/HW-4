% Parameters
fmax = 1e9;				%Max Frequency
dt = 0.5*dz/c0				%Time Step
tau = 0.5/fmax; 			%Pulse Duration


% Compute Gaussian Source
nzc = round(Nz/2)			%position of source
dt = 0.6*dz/c0
ta = [0:STEPS-1]*dt;			%time axis
t0 = 6*tau				%Pulse Position
s = dz/(2*c0) + dt/2			%Total delay between E and H
Esrc = exp(-((ta-t0/tau).^2);		%E field source
A = - sqrt(ER(nzc)/UR(nzc));		%amplitude of H field
Hsrc = A * exp(-((ta-t0+s)/tau).^2);	% H field source

