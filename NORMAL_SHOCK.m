% NORMAL SHOCK CALCULATOR
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
%   Compute the normal shock relations for a given input
%   Output solution as a structure if no output variable is defined
%   Output solution as a value if output variable is defined
%
% INPUTS
%   - inVal    : Input value [#]
%   - inVar    : Input variable [str]
%                   Upstream Mach number                = 'M1'
%                   Downstream Mach number              = 'M2'
%                   Static pressure ratio               = 'P2P1'
%                   Static temperature ratio            = 'T2T1'
%                   Static density ratio                = 'r2r1'
%                   Stagnation pressure ratio           = 'P02P01'
%                   Static-to-stagnation pressure ratio = 'P1P02'
%   - g        : Ratio of specific heats []
%   - [outVar] : Output variable desired
%                   If omitted, solution will be structure with all
%                   the output variables included
% 
% OUTPUTS
%   - sol : If [outVar] = 0, structure with the following contents
%               Upstream Mach number                = M1
%               Downstream Mach number              = M2
%               Static pressure ratio               = P2P1
%               Stagnation pressure ratio           = P02P01
%               Static density ratio                = r2r1
%               Static temperature ratio            = T2T1
%               Static-to-stagnation pressure ratio = P1P02
%            If [outVar] = 'str', will contain single value of that string
% 
% USAGE
%   Example 1: Will return a structure of all solution variables
%       sol = NORMAL_SHOCK(2.4,'M1',1.4);
%   Example 2: Will return the value of the downstream Mach number
%       sol = NORMAL_SHOCK(2.4,'M1',1.4,'M2');

function [sol] = NORMAL_SHOCK(inVal,inVar,g,outVar)

% Check input argument number
if (nargin == 3)
    outVar = 0;
end

% Catch errors associated with input variable
inVarArray = {'M1';'M2';'P2P1';'T2T1';'r2r1';'P02P01';'P1P02'};
if (~ismember(inVar,inVarArray))
    sol = inf;
    fprintf('Input variable name is incorrect!\n');
    return;
end

% User input variables
v = inVal;
if (strcmpi(inVar,'M1'))
    i = 1;
elseif (strcmpi(inVar,'M2'))
    i = 2;
elseif (strcmpi(inVar,'P2P1'))
    i = 3;
elseif (strcmpi(inVar,'T2T1'))
    i = 4;
elseif (strcmpi(inVar,'r2r1'))
    i = 5;
elseif (strcmpi(inVar,'P02P01'))
    i = 6;
elseif (strcmpi(inVar,'P1P02'))
    i = 7;
end

% Convenient parameters
gm1   = g-1;
gp1   = g+1;
gogm1 = g/gm1;
gogp1 = g/gp1;
oogm1 = 1/gm1;

% Check that specific heat ratio is greater than unity
if (g <= 1)
    sol = inf;
    fprintf('Gamma must be greater than 1\n');
    return;
end

% Solve using: M1
if (i == 1)
    if (v <= 1)
        sol = inf;
        fprintf('M1 must be greater than 1\n');
        return;
    else
        M1 = v;
    end
end

% Solve using: M2
if (i == 2)
    checkVal = sqrt((gm1)/2/g);
    if (v >= 1 || v <= checkVal)
        sol = inf;
        fprintf('M2 must be between %3.2f and 1\n',checkVal);
        return;
    else
        M1 = M2(g,v);
    end
end

% Solve using: P2/P1
if (i == 3)
    if (v <= 1)
        sol = inf;
        fprintf('P2/P1 must be greater than 1\n');
        return;
    else
        M1 = sqrt((v-1)*gp1/2/g+1);
    end
end

% Solve using: T2/T1
if (i == 4)
    if (v <= 1)
        sol = inf;
        fprintf('T2/T1 must be greater than 1\n');
        return;
    else
        aa = 2*g*gm1;
        bb = 4*g-gm1*gm1-v*gp1*gp1;
        cc = -2*gm1;
        M1 = sqrt((-bb+sqrt(bb^2-4*aa*cc))/2/aa);
    end
end

% Solve using: rho2/rho1
if (i == 5)
    checkVal = gp1/gm1;
    if (v <= 1 || v >= checkVal)
        sol = inf;
        fprintf('rho2/rho1 must be between 1 and %3.2f\n',checkVal);
        return;
    else
        M1 = sqrt(2*v/(gp1-v*gm1));
    end
end

% Solve using: P02/P01
if (i == 6)
    if (v >= 1 || v <= 0)
        sol = inf;
        fprintf('P02/P01 must be between 0 and 1\n');
        return;
    else
        Mnew = 2;
        M1   = 0;
        while (abs(Mnew-M1) > 0.00001)
            M1     = Mnew;
            al     = gp1*M1^2/(gm1*M1^2+2);
            be     = gp1/(2*g*M1^2-gm1);
            daldm1 = (2/M1-2*M1*gm1/(gm1*M1^2+2))*al;
            dbedm1 = -4*g*M1*be/(2*g*M1^2-gm1);
            fm     = (al^gogm1) * (be^oogm1) - v;
            fdm    = gogm1 * (al^oogm1) * ...
                        daldm1*(be^oogm1) + ...
                        (al^gogm1)/gm1 * ...
                        (be^((2-g)/gm1))*dbedm1;
            Mnew   = M1-(fm/fdm);
        end
    end
end

% Solve using: P1/P02
if (i == 7)
    vmax = ((g+1)/2)^(-g/(g-1));
    if (v >= vmax || v <= 0)
        sol = inf;
        fprintf('P1/P02 must be between 0 and %3.2f\n',vmax);
        return;
    else
        Mnew = 2;
        M1   = 0;
        while (abs(Mnew-M1) > 0.00001)
            M1     = Mnew;
            al     = (g+1)*M1^2/2.;
            be     = (g+1)/(2*g*M1^2-(g-1));
            daldm1 = M1*(g+1);
            dbedm1 = -4*g*M1*be/(2*g*M1^2-(g-1));
            fm     = (al^gogm1)*(be^oogm1)-(1/v);
            fdm    = gogm1 * (al^oogm1) * ...
                        daldm1*(be^oogm1) + ...
                        (al^gogm1)/gm1 * ...
                        (be^((2-g)/gm1))*dbedm1;
            Mnew   = M1-(fm/fdm);
        end
    end
end

% Solve for output variables
Mach2  = M2(g,M1);
P2P1   = 1+(2*gogp1*(M1^2-1));
P02P01 = PP0(g,M1)/PP0(g,M2(g,M1))*P2P1;
r2r1   = RR0(g,M2(g,M1))/RR0(g,M1)*P02P01;
T2T1   = TT0(g,M2(g,M1))/TT0(g,M1);
P1P02  = PP0(g,M1)/P02P01;

if (outVar == 0)
    sol.M1     = M1;
    sol.M2     = Mach2;
    sol.P2P1   = P2P1;
    sol.P02P01 = P02P01;
    sol.r2r1   = r2r1;
    sol.T2T1   = T2T1;
    sol.P1P02  = P1P02;
elseif (strcmpi(outVar,'M1'))
    sol = M1;
elseif (strcmpi(outVar,'M2'))
    sol = Mach2;
elseif (strcmpi(outVar,'P2P1'))
    sol = P2P1;
elseif (strcmpi(outVar,'P02P01'))
    sol = P02P01;
elseif (strcmpi(outVar,'r2r1'))
    sol = r2r1;
elseif (strcmpi(outVar,'T2T1'))
    sol = T2T1;
elseif (strcmpi(outVar,'P1P02'))
    sol = P1P02;
end
    

% =====================
% ===== FUNCTIONS =====
% =====================

function [M2_Out] = M2(g,M1)
    M2_Out = sqrt((1+0.5*(g-1)*M1^2)/(g*M1^2-0.5*(g-1)));
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


end



