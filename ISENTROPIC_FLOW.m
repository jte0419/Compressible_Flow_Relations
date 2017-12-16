% ISENTROPIC FLOW CALCULATOR
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
%   Compute the isentropic flow relations
%   Output solution as a structure if no output variable is defined
%   Output solution as a value if output variable is defined
% 
% INPUTS
%   - inVal    : Input value [#]
%   - inVar    : Input variable [str]
%                   Mach number                      = 'M'
%                   Static-to-stagnation temperature = 'TT0'
%                   Static-to-stagnation pressure    = 'PP0'
%                   Static-to-stagnation density     = 'rr0'
%                   Subsonic area ratio              = 'Asub'
%                   Supersonic area ratio            = 'Asup'
%                   Mach wave angle                  = 'mu'
%                   Prandtl-Meyer angle              = 'nu'
%   - g        : Ratio of specific heats []
%   - [outVar] : Output variable desired
%                   If omitted, solution will be structure with all
%                   the output variables included
% 
% OUTPUTS
%   - sol : Structure of solution variables for the given input
%               Mach number                      = M
%               Static-to-stagnation temperature = TT0
%               Static-to-stagnation pressure    = PP0
%               Static-to-stagnation density     = rr0
%               Static-to-star temperature       = TTs
%               Static-to-star pressure          = PPs
%               Static-to-star density           = rrs
%               Area ratio                       = AAs
%               Mach wave angle                  = mu
%               Prandtl-Meyer angle              = nu
%            If [outVar] = 'str', will contain single value of that string
% 
% USAGE
%   Example 1: Will return a structure of all solution variables
%       sol = ISENTROPIC_FLOW(2.6,'M',1.4);
%   Example 2: Will return the value of the area ratio
%       sol = ISENTROPIC_FLOW(2.6,'M',1.4,'AAs');

function [sol] = ISENTROPIC_FLOW(inVal,inVar,g,outVar)

% Check input argument number
if (nargin == 3)
    outVar = 0;
end

% Catch errors associated with input variable
inVarArray = {'M';'TT0';'PP0';'rr0';'Asub';'Asup';'mu';'nu'};
if (~ismember(inVar,inVarArray))
    sol = inf;
    fprintf('Input variable name is incorrect!\n');
    return;
end

% User input variables
v = inVal;
if (strcmpi(inVar,'M'))
    i = 1;
elseif (strcmpi(inVar,'TT0'))
    i = 2;
elseif (strcmpi(inVar,'PP0'))
    i = 3;
elseif (strcmpi(inVar,'rr0'))
    i = 4;
elseif (strcmpi(inVar,'Asub'))
    i = 5;
elseif (strcmpi(inVar,'Asup'))
    i = 6;
elseif (strcmpi(inVar,'mu'))
    i = 7;
elseif (strcmpi(inVar,'nu'))
    i = 8;
end

% Convenient parameters
gm1   = g-1;
gm1og = gm1/g;

% Check that specific heat ratio is greater than unity
if (g <= 1)
    fprintf('Gamma must be greater than 1\n');
    return;
end

% Solve using: Mach Number
if (i == 1)
    if (v <= 0)
        sol = inf;
        fprintf('M must be greater than 0\n');
        return;
    else
        M = v;
    end
end

% Solve using: T/T0
if (i == 2)
    if (v >= 1 || v <= 0)
        sol = inf;
        fprintf('T/T0 must be between 0 and 1\n');
        return;
    else
        M = sqrt(2*((1/v)-1)/(g-1));
    end
end

% Solve using: P/P0
if (i == 3)
    if (v >= 1 || v <= 0)
        sol = inf;
        fprintf('P/P0 must be between 0 and 1\n');
        return;
    else
        M = sqrt(2*((1/(v^gm1og))-1)/gm1);
    end
end

% Solve using: rho/rho0
if (i == 4)
    if (v >= 1 || v <= 0)
        sol = inf;
        fprintf('rho/rho0 must be between 0 and 1\n');
        return;
    else
        M = sqrt(2*((1/(v^gm1))-1)/gm1);
    end
end

% Solve using: A/A* (sub and sup)
if (i == 5 || i == 6)
    if (v <= 1)
        sol = inf;
        fprintf('A/A* must be greater than 1\n');
        return;
    else
        Mnew = 0.00001;
        M    = 0;
        if (i == 6)
            Mnew = 2;
        end
        
        while (abs(Mnew-M) > 0.000001)
            M    = Mnew;
            phi  = AAS(g,M);
            s    = (3-g)/(g+1);
            Mnew = M-(phi-v)/((phi*M)^s-phi/M);
        end
    end
end

% Solve using: Mach Angle (deg)
if (i == 7)
    if (v <= 0 || v >= 90)
        sol = inf;
        fprintf('Mach angle must be between 0 and 90 degrees\n');
        return;
    else
        M = 1/(sind(v));
    end
end

% Solve using: P-M Angle (deg)
if (i == 8)
    numax = (sqrt((g+1)/(g-1))-1)*90;
    if (v <= 0 || v >= numax)
        sol = inf;
        fprintf('P-M angle must be between 0 and %3.2f degrees\n',numax);
        return;
    else
        Mnew = 2;
        M    = 0;
        while(abs(Mnew-M) > 0.00001)
            M    = Mnew;
            fm   = (NU(g,M)-v)*(pi/180);
            fdm  = sqrt((M^2)-1)/(1+0.5*(g-1)*(M^2))/M;
            Mnew = M - (fm/fdm);
        end
    end
end

% Solve for Mach wave angle and PM angle
if (M > 1)
    mu = asind(1/M);
    nu = NU(g,M);
elseif (M == 1)
    mu = 90;
    nu = 0;
else
    mu = inf;
    nu = inf;
end

% Set solution variables
if (outVar == 0)
    sol.mu  = mu;
    sol.nu  = nu;
    sol.M   = M;
    sol.TT0 = TT0(g,M);
    sol.PP0 = PP0(g,M);
    sol.rr0 = RR0(g,M);
    sol.TTs = TTS(g,M);
    sol.PPs = PPS(g,M);
    sol.rrs = RRS(g,M);
    sol.AAs = AAS(g,M);
elseif (strcmpi(outVar,'mu'))
    sol = mu;
elseif (strcmpi(outVar,'nu'))
    sol = nu;
elseif (strcmpi(outVar,'M'))
    sol = M;
elseif (strcmpi(outVar,'TT0'))
    sol = TT0(g,M);
elseif (strcmpi(outVar,'PP0'))
    sol = PP0(g,M);
elseif (strcmpi(outVar,'rr0'))
    sol = RR0(g,M);
elseif (strcmpi(outVar,'TTS'))
    sol = TTS(g,M);
elseif (strcmpi(outVar,'PPS'))
    sol = PPS(g,M);
elseif (strcmpi(outVar,'rrs'))
    sol = RRS(g,M);
elseif (strcmpi(outVar,'AAs'))
    sol = AAS(g,M);
end

% =====================
% ===== FUNCTIONS =====
% =====================

function [nu_Out] = NU(g,M)
    term1  = sqrt((g+1)/(g-1));
    term2  = atand(sqrt(((g-1)/(g+1))*((M^2)-1)));
    term3  = atand(sqrt((M^2)-1));
    nu_Out = term1*term2 - term3;
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

function [pps_Out] = PPS(g,M)
    pps_Out = PP0(g,M)*((g+1)/2)^(g/(g-1));
end

function [rrs_Out] = RRS(g,M)
    rrs_Out = RR0(g,M)*((g+1)/2)^(1/(g-1));
end

function [tts_Out] = TTS(g,M)
    tts_Out = TT0(g,M)*((g+1)/2);
end

function [aas_Out] = AAS(g,M)
    aas_Out = (1/RRS(g,M))*sqrt(1/TTS(g,M))/M;
end

end



