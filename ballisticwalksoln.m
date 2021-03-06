function ballisticwalk
% Perform forward simulation of ballistic walking  
  
% Instructions:
%
% 1. Use the Dynamics Workbench to generate equations of motion and
%    export them to an m-file "fballwalk" for state-derivative function.
%
% 2. Copy fballwalk into the present file, as a nested function. This
%    makes the parameter values (defined below) accessible to fballwalk.
%
% 3. Write an event function "eventballwalk" to detect the knee
%    reaching full extension, and use it to terminate the ode45 integration.
%
% 4. Copy the expressions for kinetic and potential energy into a 
%    separate function energyballwalk.
%
% 5. Run a forward simulation and verify energy conservation.

% some initial states, in order q1-q3, u1-u3   
x0 = [14; -14; -60; -50; 250; -150] * pi/180; 
% above we convert deg or deg/s into rad or rad/s
   
% Set the parameter values for simulation
mshank = 0.06; mthigh = 0.097; % segment masses as fraction of body mass

lfoot = 0.25; % foot length
lshank = 0.5; lthigh = 0.5; % segment lengths, thigh and shank
lcshank = 0.437*lshank; lcthigh = 0.433*lthigh; % com locations, from
% proximal end to center of mass
lstance = lshank + lthigh;
lcstance = (mshank*(lshank-lcshank) + mthigh*(lshank+lthigh-lcthigh))/(mshank + mthigh);
% stance leg has center of mass measured from ankle

etashank = 0.735*lshank; etathigh = 0.54*lthigh; % radius of gyration about proximal end
Ishank = mshank*etashank^2 - mshank*lcshank^2; % shank inertia about com
Ithigh = mthigh*etathigh^2 - mthigh*lcthigh^2; % thigh inertia about com
Istance = Ishank + mshank*(lcstance-(lshank-lcshank))^2 + ... % shank about stance com
  Ithigh + mthigh*(-lcstance+lstance-lcthigh)^2; % thigh about stance com

%   Leg measurements
%A<-------------------------->X<------------------------>Hip
%           lcstance               lstance-lcstance

%A<------------->X<--------->K<-------------->X<-------->Hip
%  lshank-lcshank   lcshank    lthigh-lcthigh    lcthigh

g = 9.81; % gravity  
BodyMass = 65; % convert to actual mass
M1 = (mshank+mthigh)*BodyMass;
M2 = mthigh*BodyMass;
M3 = mshank*BodyMass;
I1 = Istance*BodyMass;
I2 = Ithigh*BodyMass; 
I3 = Ishank*BodyMass;
l1 = lstance; l2 = lthigh; l3 = lshank;
lc1 = lcstance; lc2 = lcthigh; lc3 = lcshank;

% Notice that the initial leg configuration does not start in double
% support.

clf; axis equal; hold on; % prepare a figure window
drawlegs(x0);
pause;


% simulate ballistic walking forward for a short duration
options = odeset('events', @eventballwalk);
[ts, xs] = ode45(@fballwalk, [0 0.5], x0);

% Here is a plot of the segment angles, defined ccw from vertical
clf; subplot(131)
plot(ts, xs(:,1:3)); 
xlabel('time'); ylabel('segment angle');
legend('q1', 'q2', 'q3');

% Now check energy conservation, as a post-processing step
energies = zeros(length(ts),1);
kineticEnergies = zeros(length(ts),1);
gravPotentialEnergies = zeros(length(ts),1);

for i = 1:length(ts)
  [energies(i), kineticEnergies(i), gravPotentialEnergies(i)] = ...
    energyballwalk(xs(i,:));
end

% show all three energies
subplot(132)
plot(ts, energies, ts, kineticEnergies, ts, gravPotentialEnergies);
xlabel('time'); ylabel('energy');
legend('total energy', 'kinetic energy', 'potential energy');

% and zoom in on just total energy
subplot(133)
plot(ts, energies);
xlabel('time'); ylabel('energy');
legend('total energy');


% Animation




% END OF MAIN FUNCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% fballwalk goes here

function xdot = fballwalk(t, x)

% State derivative code generated by Dynamics Workbench Sat 2 Feb 2019 00:22:24
% Define constants

% Define forces: 

% State assignments
q1 = x(1); q2 = x(2); q3 = x(3); 
u1 = x(4); u2 = x(5); u3 = x(6); 

s1 = sin(q1); s2 = sin(q2); s3 = sin(q3); c1m2 = cos(q1 - q2); c1m3 = cos(q1 - q3); c2m3 = cos(q2 - q3); s1m2 = sin(q1 - q2); s1m3 = sin(q1 - q3); s2m3 = sin(q2 - q3); 

% Mass Matrix
MM = zeros(3,3);
MM(1,1) = I1 + M2*(l1*l1) + M3*(l1*l1) + M1*(lc1*lc1); MM(1,2) = ...
-(c1m2*l1*lc2*M2) - c1m2*l1*l2*M3; MM(1,3) = -(c1m3*l1*lc3*M3); 
MM(2,1) = MM(1,2); MM(2,2) = I2 + M3*(l2*l2) + M2*(lc2*lc2); MM(2,3) = ...
c2m3*l2*lc3*M3; 
MM(3,1) = MM(1,3); MM(3,2) = MM(2,3); MM(3,3) = I3 + M3*(lc3*lc3); 

% righthand side terms
rhs = zeros(3,1);
rhs(1) = s1*g*lc1*M1 + s1*g*l1*M2 + s1*g*l1*M3 + s1m2*(l1*lc2*M2 + l1*lc2*M3 ...
- l1*(-l2 + lc2)*M3)*(u2*u2) + s1m3*l1*lc3*M3*(u3*u3); 
rhs(2) = -(s2*g*lc2*M2) - s2*g*l2*M3 + s1m2*(-(lc1*lc2*M2) + (-l1 + ...
lc1)*lc2*M2 - l2*(l1 - lc1)*M3 - l2*lc1*M3)*(u1*u1) - s2m3*l2*lc3*M3*(u3*u3); 
rhs(3) = -(s3*g*lc3*M3) + s1m3*(-(lc1*lc3*M3) + (-l1 + lc1)*lc3*M3)*(u1*u1) + ...
s2m3*((l2 - lc2)*lc3*M3 + lc2*lc3*M3)*(u2*u2); 

udot = MM\rhs;
xdot = [x(3+1:2*3); udot];

end % function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% energyballwalk computes kinetic and potential energy
function [energy, kineticEnergy, gravPotentialEnergy] = energyballwalk(x)
  
% State assignments
q1 = x(1); q2 = x(2); q3 = x(3); 
u1 = x(4); u2 = x(5); u3 = x(6); 

s1 = sin(q1); s2 = sin(q2); s3 = sin(q3); c1m2 = cos(q1 - q2); c1m3 = cos(q1 - q3); c2m3 = cos(q2 - q3); s1m2 = sin(q1 - q2); s1m3 = sin(q1 - q3); s2m3 = sin(q2 - q3); 

kineticEnergy = (I1*(u1*u1))/2. + (I2*(u2*u2))/2. + (I3*(u3*u3))/2. + ...
(M1*(u1*u1)*(lc1*lc1))/2. + (M2*(2*c1m2*u2*(-(u1*(l1 - lc1)) - u1*lc1)*lc2 + ...
(-(u1*(l1 - lc1)) - u1*lc1)*(-(u1*(l1 - lc1)) - u1*lc1) + ...
u2*u2*(lc2*lc2)))/2. + (M3*(2*c1m2*(-(u1*(l1 - lc1)) - u1*lc1)*(u2*lc2 - ...
u2*(-l2 + lc2)) + 2*c1m3*u3*(-(u1*(l1 - lc1)) - u1*lc1)*lc3 + ...
2*c2m3*u3*(u2*lc2 - u2*(-l2 + lc2))*lc3 + (-(u1*(l1 - lc1)) - ...
u1*lc1)*(-(u1*(l1 - lc1)) - u1*lc1) + (u2*lc2 - u2*(-l2 + lc2))*(u2*lc2 - ...
u2*(-l2 + lc2)) + u3*u3*(lc3*lc3)))/2.;

gravPotentialEnergy = g*lc1*M1*cos(q1) - M2*(-(g*l1*cos(q1)) + g*lc2*cos(q2)) ...
- M3*(-(g*l1*cos(q1)) + g*l2*cos(q2) + g*lc3*cos(q3));


energy = kineticEnergy + gravPotentialEnergy;

end % energyballwalk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value, direction, isterminal] = eventballwalk(x)

% WRITE A FUNCTION TO DETECT KNEE GETTING STRAIGHT

% State assignments
q1 = x(1); q2 = x(2); q3 = x(3); 
u1 = x(4); u2 = x(5); u3 = x(6); 

% swing knee angle (zero at full extension)
  value = q2 - q3;

  direction = -1; % from above

  isterminal = 1; % stop the integration

end % eventballwalk

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function drawlegs(x);
    % Use the angles to draw a stick figure of the leg configuration
    % The segment lengths should be defined in outer scope of this
    % nested function.
    % A figure window should already exist, with axis equal 
    % recommended.
    
    q1 = x(1); q2 = x(2); q3 = x(3); 

    xtoe = 0;
    xankle = xtoe - lfoot;
    xhip = xankle - lstance*sin(q1);
    xswingknee = xhip + lthigh*sin(q2);
    xswingankle = xswingknee + lshank*sin(q3);
    xswingtoe = xswingankle + lfoot*cos(q3);
    ytoe = 0;
    yankle = ytoe;
    yhip = yankle + lstance*cos(q1);
    yswingknee = yhip - lthigh*cos(q2);
    yswingankle = yswingknee - lshank*cos(q3);
    yswingtoe = yswingankle - lfoot*cos(q3);

    xlocations = [xtoe, xankle, xhip, xswingknee, xswingankle, xswingtoe];
    ylocations = [ytoe, yankle, yhip, yswingknee, yswingankle, yswingtoe];

    plot(xlocations, ylocations, 'linewidth', 3); % segment lines
    plot(xlocations, ylocations, '.', 'markersize', 10 ); % a dot at each joint
end % drawlegs

end % ballisticwalk 
