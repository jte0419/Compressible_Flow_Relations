% OBLIQUE SHOCK CALCULATOR
% Adapted by  : JoshTheEngineeer
% Website     : www.JoshTheEngineer.com
% YouTube     : www.youtube.com/JoshTheEngineer
% Based on    : VT Calculator
%               http://www.dept.aoe.vt.edu/~devenpor/aoe3114/calc.html
% Developed by: William J. Devenport, Virginia Tech
%               Adam Ford (Fanno Flow, Rayleigh Flow)
%               Stephen Krauss (Conical Flow)
% Started: 12/13/17
% Updated: 12/13/17 - Copied over from GUI
%                   - Works as intended
%          12/14/17 - Added input variable name check
% 
% PUPROSE
%   Compute the oblique shock relations for a given input
%   Output solution as a structure if no output variable is defined
%   Output solution as a value if output variable is defined
%
% INPUTS
%   - inVal    : Input value [#]
%   - inVar    : Input variable [str]
%                   Turn angle (weak shock)     = 'thetaWeak'
%                   Turn angle (strong shock)   = 'thetaStrong'
%                   Shock angle                 = 'beta'
%                   Upstream normal Mach number = 'Mn1'
%   - M1       : Upstream Mach number []
%   - g        : Ratio of specific heats []
%   - [outVar] : Output variable desired
%                   If omitted, solution will be structure with all
%                   the output variables included
% 
% OUTPUTS
%   - sol : If [outVar] = 0, structure with the following contents
%               Turn angle                    = theta
%               Shock angle                   = beta
%               Upstream Mach number          = M1
%               Upstream normal Mach number   = Mn1
%               Downstream normal Mach number = Mn2
%               Downstream Mach number        = M2
%               Static pressure ratio         = P2P1
%               Stagnation pressure ratio     = P02P01
%               Static density ratio          = r2r1
%               Static temperature ratio      = T2T1
%               Shock state                   = state
%            If [outVar] = 'str', will contain single value of that string
% 
% USAGE
%   Example 1: Will return a structure of all solution variables
%       sol = OBLIQUE_SHOCK(20,'thetaWeak',2.5,1.4);
%   Example 2: Will return the value of the static pressure ratio
%       sol = OBLIQUE_SHOCK(20,'thetaWeak',2.5,1.4,'P2P1');

function [sol] = OBLIQUE_SHOCK(inVal,inVar,M1,g,outVar)

% Check input argument number
if (nargin == 4)
    outVar = 0;
end

% Catch errors associated with input variable
inVarArray = {'thetaWeak';'thetaStrong';'beta';'Mn1'};
if (~ismember(inVar,inVarArray))
    sol = inf;
    fprintf('Input variable name is incorrect!\n');
    return;
end

% User input variables
v = inVal;
if (strcmpi(inVar,'thetaWeak'))
    i = 1;
elseif (strcmpi(inVar,'thetaStrong'))
    i = 2;
elseif (strcmpi(inVar,'beta'))
    i = 3;
elseif (strcmpi(inVar,'Mn1'))
    i = 4;
end

% Check that specific heat ratio is greater than unity
if (g <= 1)
    sol = inf;
    fprintf('Gamma must be greater than 1\n');
    return;
end

% Make sure Mach number is greater than unity
if (M1 <= 1)
    sol = inf;
    fprintf('M1 must be greater than 1\n');
    return;
end

% Solve using: Turn Angles (Weak and Strong)
if (i == 1 || i == 2)
    theta = v;
    if (theta >= 90)
        sol = inf;
        fprintf('Turning angle is too large\n');
        return;
    end
    
    % Check to make sure turning angle is greater than zero
    if (theta <= 0)
        sol = inf;
        fprintf('Turning angle must be greater than zero\n');
        return;
    end
    
    % Solve for shock angle [deg]
    beta = MDB(g,M1,theta,i);
	
    % If shock angle is out of bounds
    if (beta < 0)
        if (outVar == 0)
            sol.beta   = inf;
            sol.theta  = inf;
            sol.Mn1    = inf;
            sol.Mn2    = inf;
            sol.M2     = inf;
            sol.P2P1   = inf;
            sol.P02P01 = inf;
            sol.r2r1   = inf;
            sol.T2T1   = inf;
            sol.state  = 'Detached';
        else
            sol = 'Detached';
        end
        return;
    end

% Solve using: Shock Angle
elseif (i == 3)
    beta = v;
    if (beta >= 90)
        sol = inf;
        fprintf('Shock angle must be less than 90 degrees\n');
        return;
    end
    
    if (beta-asind(1/M1) <= 0)
        sol = inf;
        fprintf('Shock angle must be greater than Mach wave angle\n');
        return;
    end
    
    theta = MBD(g,M1,beta);
    
% Solve using: M1n
elseif (i == 4)
    Mn1 = v;
    
    if (Mn1 <= 1 || Mn1 >= M1)
        sol = inf;
        fprintf('Mn1 must be between 1 and M1\n');
        return;
    end
    
    beta  = asind(Mn1/M1);                                                  % Shock angle [deg]
    theta = MBD(g,M1,beta);                                                 % Turn angle [deg]
end

% Solve for solution variables
Mn1      = M1*sind(beta);                                                   % Upstream normal Mach number []
Mn2      = M2(g,Mn1);                                                       % Downstream normal Mach number []
Mach2    = Mn2/sind(beta-theta);                                            % Downstream Mach number []
P2P1     = 1+2*g/(g+1)*(Mn1^2-1);                                           % Static pressure ratio []
P02P01   = PP0(g,Mn1)/PP0(g,M2(g,Mn1))*P2P1;                                % Stagnation pressure ratio []
rho2rho1 = RR0(g,M2(g,Mn1))/RR0(g,Mn1)*P02P01;                              % Static density ratio []
T2T1     = TT0(g,M2(g,Mn1))/TT0(g,Mn1);                                     % Static temperature ratio []

% Set solution variables
if (outVar == 0)                                                            % If no output variable is defined
    sol.theta  = theta;
    sol.beta   = beta;
    sol.M1     = M1;
    sol.Mn1    = Mn1;
    sol.Mn2    = Mn2;
    sol.M2     = Mach2;
    sol.P2P1   = P2P1;
    sol.P02P01 = P02P01;
    sol.r2r1   = rho2rho1;
    sol.T2T1   = T2T1;
    sol.state  = 'Attached';
elseif (strcmpi(outVar,'theta'))                                            % Output turn angle [deg]
    sol = theta;
elseif (strcmpi(outVar,'beta'))                                             % Output shock angle [deg]
    sol = beta;
elseif (strcmpi(outVar,'M1'))                                               % Output upstream Mach number []
    sol = M1;
elseif (strcmpi(outVar,'Mn1'))                                              % Output upstream normal Mach number []
    sol = Mn1;
elseif (strcmpi(outVar,'Mn2'))                                              % Output downstream normal Mach number []
    sol = Mn2;
elseif (strcmpi(outVar,'M2'))                                               % Output downstream Mach number []
    sol = Mach2;
elseif (strcmpi(outVar,'P2P1'))                                             % Output static pressure ratio []
    sol = P2P1;
elseif (strcmpi(outVar,'P02P01'))                                           % Output stagnation pressure ratio []
    sol = P02P01;
elseif (strcmpi(outVar,'r2r1'))                                             % Output static density ratio []
    sol = rho2rho1;
elseif (strcmpi(outVar,'T2T1'))                                             % Output static temperature ratio[]
    sol = T2T1;
elseif (strcmpi(outVar,'state'))                                            % Output whether shock is attached or detached
    sol = 'Attached';
end

% =====================
% ===== FUNCTIONS =====
% =====================

function [M2_Out]  = M2(g,M1)
    M2_Out = sqrt((1+0.5*(g-1)*M1*M1)/(g*M1*M1-0.5*(g-1)));
end

function [pp0_Out] = PP0(g,M)
    pp0_Out = (1+(g-1)/2*(M^2))^(-g/(g-1));
end

function [rr0_Out] = RR0(g,M)
    rr0_Out = (1+(g-1)/2*(M^2))^(-1/(g-1));
end

function [tt0_Out] = TT0(g,M)
    tt0_Out = (1+(g-1)/2*(M^2))^(-1);
end

function [mbd_Out] = MBD(g,M1,b)
    mbd_Out = atand((M1^2*sind(2*b)-2/tand(b))/(2+M1^2*(g+cosd(2*b))));
end

function [mdb_Out] = MDB(g,M1,d,i)
    d = d*(pi/180);
    p = -(M1*M1+2)/M1/M1-g*sin(d)*sin(d);
    q = (2*M1*M1+1)/(M1^4)+((g+1)*(g+1)/4+(g-1)/M1/M1)*sin(d)*sin(d);
    r = -cos(d)*cos(d)/(M1^4);
    
    a = (3*q-p*p)/3;
    b = (2*p*p*p-9*p*q+27*r)/27;
    
    test = b*b/4+a*a*a/27;
    
    if (test > 0)
        mdb_Out = -1;
        return;
    else
        if (test == 0)
            x1 = sqrt(-a/3);
            x2 = x1;
            x3 = 2*x1;
            if (b > 0)
                x1 = -1*x1;
                x2 = -1*x2;
                x3 = -1*x3;
            end
        end

        if (test < 0)
            phi = acos(sqrt(-27*b^2/4/a/a/a));
            x1  = 2*sqrt(-a/3)*cos(phi/3);
            x2  = 2*sqrt(-a/3)*cos(phi/3+pi*2/3);
            x3  = 2*sqrt(-a/3)*cos(phi/3+pi*4/3);
            if (b > 0)
                x1 = -1*x1;
                x2 = -1*x2;
                x3 = -1*x3;
            end
        end

        s1 = x1-p/3;
        s2 = x2-p/3;
        s3 = x3-p/3;

        if (s1 < s2 && s1 < s3)
          t1 = s2;
          t2 = s3;
        elseif (s2 < s1 && s2 < s3)
          t1 = s1;
          t2 = s3;
        else
          t1 = s1;
          t2 = s2;
        end
        
        b1 = asin(sqrt(t1));
        b2 = asin(sqrt(t2));
        
        betas = b1;
        betaw = b2;
        if (b2 > b1)
          betas = b2;
          betaw = b1;
        end
        
        if (i == 1)
            mdb_Out = betaw*(180/pi);
            return;
        end
        
        if (i == 2)
            mdb_Out = betas*(180/pi);
            return;
        end
    end
end


end



