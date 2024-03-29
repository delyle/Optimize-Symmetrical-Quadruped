function output = SymQuadOptCtrl_NCW(auxdata,guess)

%-------------------------------------------------------------------%
%-------------------- Data Required by Problem ---------------------%
%-------------------------------------------------------------------%

% specify auxdata if not already done

D = auxdata.D;
l = auxdata.lmax;
if ~isfield(auxdata,'abounds')
    auxdata.abounds = [-pi/2 pi/2];
end
abounds = auxdata.abounds;

if ~isfield(auxdata,'Fmax')
    auxdata.Fmax = 10;
end
Fmax = auxdata.Fmax;

if ~isfield(auxdata,'Fdotmax')
    auxdata.Fdotmax = 100;
end
dFmax = auxdata.Fdotmax;

if ~isfield(auxdata,'meshtolerance')
    auxdata.meshtolerance = 1e-4;
end

if ~isfield(auxdata,'setMesh')
    auxdata.setMesh = [4 4];
end


%-------------------------------------------------------------------%
%------------------------- Variable Bounds -------------------------%
%-------------------------------------------------------------------%
% ----- PHASE 1 ----- %
i = 1;
bounds.phase(i).initialtime.lower = 0;              % scalar
bounds.phase(i).initialtime.upper = 0;              % scalar
bounds.phase(i).finaltime.lower = 0.5;                % scalar
bounds.phase(i).finaltime.upper = 0.5;                % scalar

% States:
% 6 kinematic states: x, y, a, u, v, w
% 5 forces: FLH FLFt FLFl FRH FRF
% 1 integrated force: int_0^t FLF dt
xlow = 0; xupp = D/2;
ylow = 0; yupp = 4*max(l);
alow = abounds(1); aupp = abounds(2);
ulow = -4*D; uupp = 4*D;
vlow = -4*D; vupp = 4*D;
wlow = -4  ; wupp = 4  ;
F_low = zeros(1,5);
Fi_upp = [0 1 0 1 1]*Fmax;
Ff_upp = [1 0 1 0 1]*Fmax;
F_upp  = [1 1 1 1 1]*Fmax;
intF_low = 0;
intF_upp = Fmax*0.5;

% Controls:
% 5 Force rate: Fdot
% 2 slack variables for Net Work: p, q
% 5 relaxation parameters for limb length constraints
% 2 relaxation parameter for simultaneous limb contact constraints
% Not worrying about torque constraints... for now!
dF_low = [1 1 1 1 1]*(-dFmax);
dF_upp = [1 1 1 1 1]*dFmax;
par_low = zeros(1,9);
par_upp = 1000*ones(1,9);

bounds.phase(i).initialstate.lower = [xlow ylow alow ulow vlow wlow F_low intF_low];           % row vector, length = numstates
bounds.phase(i).initialstate.upper = [xlow yupp aupp uupp vupp wupp Fi_upp intF_low];           % row vector, length = numstates
bounds.phase(i).state.lower = [xlow ylow alow ulow vlow wlow F_low intF_low];                  % row vector, length = numstates
bounds.phase(i).state.upper = [xupp yupp aupp uupp vupp wupp F_upp intF_upp];                  % row vector, length = numstates
bounds.phase(i).finalstate.lower = [xupp ylow alow ulow vlow wlow F_low intF_low];             % row vector, length = numstates
bounds.phase(i).finalstate.upper = [xupp yupp aupp uupp vupp wupp Ff_upp intF_upp];             % row vector, length = numstates

bounds.phase(i).control.lower = [dF_low par_low];                % row vector, length = numstates
bounds.phase(i).control.upper = [dF_upp par_upp];                % row vector, length = numstates

% 1 integral for work
% 1 integral for force rate
% 1 integral for slack penalties
bounds.phase(i).integral.lower = 0*ones(1,3);                 % row vector, length = numintegrals
bounds.phase(i).integral.upper = 100*ones(1,3);                 % row vector, length = numintegrals

% Footfall locations: PLH PTF
bounds.parameter.lower = -(1+sum(l))*[1 1];                      % row vector, length = numintegrals
bounds.parameter.upper = (D+1+sum(l))*[1 1];                      % row vector, length = numintegrals

% Endpoint constraints
% 5 Kinematic periodicity: yauvw(0.5) - yauvw(0) = 0
% 3 Force continuity: FLH(0.5) - FRH(0) = FRF(0.5) - FTF(0) = FLF(0.5) - FRF(0) = 0
%

bounds.eventgroup.lower = zeros(1,8); % row vector
bounds.eventgroup.upper = zeros(1,8); % row vector

% Path constraints
% 1 complementary: FTF*int_0^t FLF dt = 0
% 1 complementary: FLFt*FLFldot = 0
% 5 complementary: F_i*(l_i - l_imax) >= 0
% 2 normal: y - rH*sin(a) >= 0, y + rF*sin(a) >= 0
% 1 normal: dot(F_net,v) - p + q = 0
% ----- PHASE 1 ----- %
i = 1;
bounds.phase(i).path.lower = zeros(1,10); % row vector, length = number of path constraints in phase
bounds.phase(i).path.upper = [0, 0, Fmax*max(l)*[1 1 1 1 1], 1 + 4*max(l)*[1 1], 0]; % row vector, length = number of path constraints in phase
%-------------------------------------------------------------------------%
%---------------------------- Provide Guess ------------------------------%
%-------------------------------------------------------------------------%
% ----- PHASE 1 ----- %

if isempty(guess)
    i = 1;
    ntime = 16;
    guess.phase(i).time    = linspace(0,0.5,ntime)';                % column vector, min length = 2
    guess.phase(i).state   = bounds2randguess(bounds.phase.state,ntime);                % array, min numrows = 2, numcols = numstates
    guess.phase(i).control = bounds2randguess(bounds.phase.control,ntime);                % array, min numrows = 2, numcols = numcontrols
    guess.phase(i).integral = bounds2randguess(bounds.phase.integral,1);               % scalar
    
    guess.parameter = bounds2randguess(bounds.parameter,1);                    % row vector, numrows = numparams
elseif isstruct(guess)
    % it's an output struct from a previous trial
    guess1 = guess;
    guess = [];
    % pull out the guess
    guess.phase.time = guess1.result.solution.phase.time;
    guess.phase.state = guess1.result.solution.phase.state;
    guess.parameter = guess1.result.solution.parameter;
    guess.phase.control = guess1.result.solution.phase.control;
    guess.phase.integral = guess1.result.solution.phase.integral;
    if isfield(auxdata,'downsample') && ~auxdata.downsample
        % Normally GPOPS will automatically downsample the guess.
        % This line ensures that is uses the exact same mesh of the guess
        setup.mesh.phase = guess1.result.setup.mesh.phase;
        auxdata.setMesh = [];
    end
    if isfield(auxdata,'perturb') && auxdata.perturb
       f = auxdata.perturb;
       guess.phase.state = addnoise(guess.phase.state,f);
       guess.parameter = addnoise(guess.parameter,f);
       guess.phase.control = addnoise(guess.phase.control,f);
       guess.phase.integral = addnoise(guess.phase.integral,f);
       guess1.result.solution = guess;
       All_results_plotSYM(guess1);
       pause(2)
    end
    clearvars guess1
end

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
if isfield(auxdata,'mesh')
   setup.mesh = auxdata.mesh; 
end

setup.mesh.maxiterations                = auxdata.meshMaxIter;
setup.mesh.tolerance                    = auxdata.meshtolerance;

if ~isempty(auxdata.setMesh)
    setup.mesh.phase.fraction = (1/auxdata.setMesh(1))*ones(1,auxdata.setMesh(1));
    setup.mesh.phase.colpoints = auxdata.setMesh(2)*ones(1,auxdata.setMesh(1));
end
%-------------------------------------------------------------------%
%--------------------------- Problem Setup -------------------------%
%-------------------------------------------------------------------%
if ~isfield(auxdata,'name')
    auxdata.name = 'SymQuadOptCtrl';
end
setup.name                        = auxdata.name;
setup.functions.continuous        = @Continuous;
setup.functions.endpoint          = @Endpoint;
setup.auxdata                     = auxdata;
setup.bounds                      = bounds;
setup.guess                       = guess;


setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'first';
setup.method = 'RPM-integration';
setup.scales.method = 'automatic-hybridupdate';

% NLP
setup.nlp.solver = 'snopt';
setup.nlp.snoptoptions.maxiterations = auxdata.snoptIter;
setup.nlp.snoptoptions.tolerance = auxdata.snopttolerance;


%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);
end


function phaseout = Continuous(input)

% extract data
%t = input.phase(1).time;
X = input.phase(1).state;
U = input.phase(1).control;
P = input.phase(1).parameter;
aux = input.auxdata;
D = aux.D;
Fr = aux.Fr;
I = aux.I;
lmax = aux.lmax;

x = X(:,1);
y = X(:,2);
a = X(:,3);
u = X(:,4);
v = X(:,5);
w = X(:,6);

F = X(:,7:11);
%FLH = F(:,1);
FLFt = F(:,2);
FLFl = F(:,3);
%FRH = F(:,4);
%FRF = F(:,5);

dF = U(:,1:5);
dFLFl = dF(:,3);
p = U(:,6);
q = U(:,7);
slimb = U(:,8:12);
ssimult = U(:,13:14);

intFLFl = X(:,12);
mf = input.auxdata.mf;
mh = (1-mf);
lb_hat =[cos(a), sin(a)];
ntime = length(x);
zs = zeros(ntime,1);
os = ones(ntime,1);
rF = mh*lb_hat;
rH = -mf*lb_hat;
xy = [x,y];
xyF = xy + rF;
xyH = xy + rH;
P1vec = [P(:,1),zs];
P2vec = [P(:,2),zs];
Dvec = [D*os,zs];
[L,Fvec] = deal(zeros(ntime,3,5));
L(:,1:2,1) = xyH - P1vec; % LLH
L(:,1:2,2) = xyF - P2vec; % LFT
L(:,1:2,3) = xyF - (P2vec + Dvec); % LFL
L(:,1:2,4) = xyH - (P1vec - Dvec/2);
L(:,1:2,5) = xyF - (P2vec + Dvec/2);

absL = sqrt(dot(L,L,2));

F3 = permute(F,[1 3 2]);

Fvec(:,1:2,:) = [F3,F3].*L(:,1:2,:)./[absL,absL];
FHtot = sum(Fvec(:,:,[1 4]),3);
FFtot = sum(Fvec(:,:,[2 3 5]),3);

Atrans = (FHtot + FFtot - [zs,os,zs])/Fr;
TorqueF = cross([rF,zs],FFtot);
TorqueH = cross([rH,zs],FHtot);
Arot = (TorqueF + TorqueH)/(I*Fr);


xdot = [u,v,w,Atrans+Arot,dF,FLFl]; % provide derivative
phaseout.dynamics = xdot;


% Path constraints
% 1 complementary: FTF*int_0^t FLF dt = 0
% 5 complementary: F_i*(l_imax-l_i) >= 0
% 2 normal: y - rH*sin(a) >= 0, y + rF*sin(a) >= 0
% 5 normal: F_i*dl - p_i + q_i = 0

SimultContact(:,1) = FLFt.*intFLFl - ssimult(:,1);
SimultContact(:,2) = FLFt.*dFLFl - ssimult(:,2);

absL = permute(absL,[1 3 2]);
LengthDiff = absL;
LengthDiff(:,[1 4]) = lmax(2) - absL(:,[1 4]); % hind
LengthDiff(:,[2 3 5]) = lmax(1) - absL(:,[2 3 5]) ; % fore
LimbLength = F.*LengthDiff - slimb;

AboveGround = [y + rH(:,2),y + rF(:,2)];

uv = [u,v];
Ftot = sum(Fvec,3);
Abs2Slack = dot(Ftot(:,1:2),uv,2) - p + q;

phaseout.path = [SimultContact,LimbLength,AboveGround,Abs2Slack]; % path constraints, matrix of size num collocation points X num path constraints


% Cost function integrand
c = aux.c;
c_Fdot = c(1);
c_pos = c(2);
c_neg = c(3);
c_pq = c(4);
c_simult = c(5);
c_limb = c(6);
phaseout.integrand = [c_pos*p + c_neg*q,...
    c_Fdot.*sum(dF.^2,2),...
    c_pq.*p.*q+c_simult.*sum(ssimult,2) + c_limb.*sum(slimb,2)];
end

function output = Endpoint(input)
% Endpoint constraints
% 5 Kinematic periodicity: yauvw(0.5) - yauvw(0) = 0
% 3 Force continuity: FLH(0.5) - FRH(0) = FRF(0.5) - FTF(0) = FLF(0.5) - FRF(0) = 0
Xi = input.phase.initialstate;
Xf = input.phase.finalstate;

Fi = Xi(7:11);
Ff = Xf(7:11);

KinPeriod = Xf(2:6) - Xi(2:6);
ForceCont = Ff([1 5 3]) - Fi([4 2 5]);


output.eventgroup.event = [KinPeriod,ForceCont]; % event constraints (row vector)
output.objective = sum(input.phase.integral); % objective function (scalar)
end

function guess = bounds2randguess(bounds,n)
rng('shuffle')
os = ones(n,1);
nbounds = length(bounds.lower);
guess = os*bounds.lower + rand(n,nbounds).*(os*(bounds.upper-bounds.lower));

end

function A = addnoise(a,f)
rng('shuffle')

s = range(a)*f;

n = size(a);

G = 1+(ones(n(1),1)*s).*randn(n);

A = G.*a;

end