function qs = assemblebwroot
% Finds the initial conditions for leg angles q1, q2, q3 that produce
% the desired step length and forward progression of the hip.
% Desired step length is 0.68 m, and forward progression of the hip
% refers to the amount the hip advances forward during double support,
% set to an amount 0.9*(foot length). Both of these values are set from
% human data, and the corresponding leg angles are solved for numerically.

% The q's are defined as q1: stance leg, q2: swing thigh, q3: swing
% shank. These are defined as positive in the counter-clockwise
% direction, measured from the vertical.

% The problem being solved is that the initial conditions do not 
% necessarily satisfy kinematic constraints, such as ensuring that
% the ballistic walking model starts with the trailing leg toe on
% the ground. The solution is numerically solved for a combination of
% the q's that satisfies this constraint.

% This version performs root finding to determine a unique assembly.
% The stance leg angle is set as an explicit constraint, so that
% the remaining two angles q2 and q3 are uniquely determined by the
% horizontal and vertical positions of the swing toe.

% YOU WILL NEED TO ADD CODE IN THREE LOCATIONS
% The rest of the program demonstrates a few other kinematic
% issues, but require no coding.

% Following parameter values are available in all subfunctions:
%  SL SPEED FWDP  lfoot l1 l2 l3

lfoot = 0.25; % foot length
lshank = 0.5; lthigh = 0.5; % segment lengths, thigh and shank
l1 = lshank + lthigh; L = l1; lstance = l1;  % stance leg
l2 = lthigh;
l3 = lshank; 
FWDP = 0.9*lfoot;   % forward progression of hip by end of double support

% The following constraints correspond to typical human walking
SL = 0.68;          % step length in m
SPEED = 1.25;       % nominal speed of 1.25 m/s

% This is the equivalent of Mochon & McMahon's equation (8):
%   SL/2 - L*sin(q10) = FWDP
q10 = asin((SL/2 - FWDP)/L);

q20 = -30*pi/180;
q30 = -50*pi/180; % initial guesses, converted from deg
q23 = [q20; q30]; % contains values for q2 and q3

% Since q1 is already determined, the error in toe position only depends
% on q2 and q3. The error returned has two components, x and y
initialToeError = toeerror(q23); 
fprintf(1,'Initial guess q = [%g %g %g]\n', [q10 q20 q30]);
fprintf(1,'Initial toe error = %g %g\n', initialToeError);

% Plot the initial and target toe positions along with a stick figure
clf; axis equal; hold on; % prepare a figure window
plot(-initialToeError(1)-SL+0*lfoot,-initialToeError(2), 'b.', ...
  -SL+0*lfoot, 0, 'r.'); hold on; % target indicated by red dot
legend('Initial toe position', 'Target toe position');
drawballwalk([q10 q20 q30]); % draw a stick model 

% now solve for the q2 and q3 that yield the desired toe position
q23star = findroot(@toeerror, q23); %% PUT YOUR CODE HERE TO USE FINDROOT TO DETERMINE ASSEMBLY 

qs = [q10; q23star(1); q23star(2)]; % return the solution

% Plot the final configuration after finding the root
fprintf(1,'Final value   q = [%g %g %g]\n', [q10 q23star(1) q23star(2)]);
fprintf(1,'Final toe error = %g %g\n', toeerror(q23star));
drawballwalk(qs);

% next part: find combinations of q's to produce desired x's

% Here is the partial derivative of the toe location with respect
% to the q's, evaluated at the assembled configuration found above.
% (Note that we can find the gradient of the toe error as equivalent
% to the gradient of the actual toe location, since the target location
% is constant. However, this uses the negative of the actual toe location
% so we will have to use -J later.
J = fjacobian(@toeerror, qs);

[U, S, V] = svd(J); % singular value decomposition

qsoln = -J \ [1; 0]; % solve for a purely horizontal motion to the right

J2 = J(2,:);        % second row of J, representing vertical motion
[U2, S2, V2] = svd(J2);

% The nullspace of the second row represents combinations of qs
% that produce zero vertical motion
nullspaceJ2 = [V2(:,2) V2(:,3)]; 

fprintf(1,'A formula for all motions that move the toe in the [1,0] velocity\n');
fprintf(1,'is  [ q1dot ]   [ %+5.3f ]        [ %+5.3f ]\n', qsoln(1), V(1,3) );
fprintf(1,'    [ q2dot ] = [ %+5.3f ] + c1 * [ %+5.3f ]\n', qsoln(2), V(2,3) );
fprintf(1,'    [ q3dot ]   [ %+5.3f ]        [ %+5.3f ]\n', qsoln(3), V(3,3) );

fprintf(1,'A formula for all motions that have no vertical toe motion is\n');
fprintf(1,'[ dq1 ]        [ %+5.3f ]        [ %+5.3f ]\n', V2(1,2), V2(1,3));
fprintf(1,'[ dq2 ] = c1 * [ %+5.3f ] + c2 * [ %+5.3f ]\n', V2(2,2), V2(2,3));
fprintf(1,'[ dq3 ]        [ %+5.3f ]        [ %+5.3f ]\n', V2(3,2), V2(3,3)); 

% Another way to do the same iteration as in the jacobian method
% is to use our particular findroot function, which performs a Newton 
% search using the pseudo-inverse of the gradient. This yields a 
% configuration that minimizes the sum-squared deviations in the q's
% with each iteration. (Note that in other methods of root finding,
% it is not necessarily possible to find the root of a function with
% two elements, if the input has three elements as is the case here.)
qsoln2 = findroot(@toeerror, [q10; q20; q30]);
fprintf(1,'The solution from findroot is [%g; %g; %g]\n', qsoln2);

%% Subfunctions below are called by main program, and have access
% to parameter values

function err = toeerror(q23) 
% Returns the horizontal and vertical toe positions p = [xtoe; ytoe],
% as an error relative to the desired toe position. The desired toe 
% position puts the toe one step length (minus a foot
% length behind the stance heel, at a height of 0 above the ground.
% If both feet were flat on the ground, the distance between heels 
% would be exactly one step length.
% We assume the foot is 90 degrees from the shank. 
% The input vector x contains q2 and q3 angles, and q1 is derived
% from a kinematic constraint on where the hip is at the end of double
% support


% This is the equivalent of Mochon & McMahon's equation (8):
%   SL/2 - L*sin(q10) = FWDP
q10 = asin((SL/2 - FWDP)/L);

q20 = q23(1); q30 = q23(2);

% This is the position of the swing toe relative to the stance heel
xtoe = (-l1*sin(q10) + l2*sin(q20) + l3*sin(q30) + lfoot*sin(q30+pi/2));
ytoe = l1*cos(q10) - l2*cos(q20) - l3*cos(q30) - lfoot*cos(q30+pi/2);

% Calculate the error between actual and desired toe positions
err = [ (-SL+lfoot)-xtoe ; 0-ytoe ]; %% PUT YOUR CODE HERE TO CALCULATE THE ERROR IN X AND Y POSITIONS:
% (This should be a vector that, when placed at the current toe position,
% will point towards the desired toe position.)

end


function hlegout = drawballwalk(x, hlegs);
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
    yswingtoe = yswingankle + lfoot*sin(q3);

    xlocations = [xtoe, xankle, xhip, xswingknee, xswingankle, xswingtoe];
    ylocations = [ytoe, yankle, yhip, yswingknee, yswingankle, yswingtoe];

    if nargin < 2 % need to produce new line
        hlegout = line('xdata', xlocations, 'ydata', ylocations, 'linewidth', 3, 'linestyle', '-',...
            'marker', '.', 'markersize', 20);
    else          % lines exist, so update their positions
        set(hlegs,'xdata', xlocations, 'ydata', ylocations);
    end
        
end % drawballwalk


function [xstar, cnvrg] = findroot(f, x0, parms)
% FINDROOT  Finds the root of a vector function, with vector-valued x0.
%  
% xstar = findroot(f, x0) performs a Newton search and returns xstar,
%   the root of the function f, starting with initial guess x0. The
%   function will typically be expressed with a function handle, e.g. 
%   @f.
%
%   Optional input findroot(f, x0, parms) includes parameters structure,
%   with fields dxtol (smallest allowable dx) and maxiter (max #
%   iterations).
%
%   Optional second output [xstar, cnvrg] = findroot... signals
%   successful convergence. It is set to 
%   false if maximum iterations is exceeded.

if nargin < 3
    parms.dxtol = 1e-6;      % default tolerance for min change in x
    parms.dftol = 1e-6;      % default tolerance for min change in f
    parms.maxiter = 1000;    % allow max of 1000 iterations before quitting
    parms.finitediffdx = []; % finite differencing step size
end

% We will take steps to improve on x, while monitoring
% the change dx, and stopping when it becomes small
dx = Inf;          % initial dx is large
iter = 0;          % count the iterations we go through
x = x0(:);         % start at this initial guess (treat x as a column vector)
fprevious = f(x0); % use to compare changes in function value
df = Inf;          % change in f is initialized as large

% Loop through the refinements, checking to make sure that x and f
% change by some minimal amount, and that we haven't exceeded
% the maximum number of iterations

while max(abs(dx)) >= parms.dxtol & max(abs(df)) >= parms.dftol & ...
        iter <= parms.maxiter
  % Here is the main Newton step
  J = fjacobian(f, x, parms.finitediffdx); % This is the slope
  dx = J\(0 - f(x));                 % and this is the correction to x
  x = x + dx;

  % Update information about changes
  iter = iter + 1;
  df = f(x) - fprevious;
  fprevious = f(x);
end

xstar = x; % use the latest, best guess

cnvrg = true;

if iter > parms.maxiter % we probably didn't find a good solution
  warning('Maximum iterations exceeded in findroot');
  cnvrg = false;
end

if size(x0,2) > 1   % x0 was given to us as a row vector
    xstar = xstar'; % so return xstar in the same shape
end
    
end % findroot

function dfdx = fjacobian(f, x0, dx)
% FJACOBIAN   Computes Jacobian (partial derivative) of function
%
% dfdx = fjacobian(@f, x0 [, dx])
%
%   Uses finite differences to compute partial derivative of vector 
%   function f evaluated at vector x0. 
%   Optional argument dx specifies finite difference
%   step. Note that argument f should typically be entered as a
%   function handle, e.g. @f.
%
%   When dx is not given or empty, a default of 1e-6 is used.

  if nargin < 3 || isempty(dx)
    dx = 1e-6;
  end
  
  f0 = f(x0);
  J = zeros(length(f0), length(x0));
  for i = 1:length(x0)
    xperturbed = x0;
    xperturbed(i) = xperturbed(i) + dx;
    df(:,i) = f(xperturbed) - f0;
  end
  dfdx = df / dx;
end

end % assembleBWroot