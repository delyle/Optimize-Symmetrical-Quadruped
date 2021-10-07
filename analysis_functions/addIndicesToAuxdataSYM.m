function output = addIndicesToAuxdataSYM(output)
% NOTE: this function only works for the axial limb work configuration!

auxdata = output.result.setup.auxdata;


% States:
% 6 kinematic states: x, y, a, u, v, w
% 5 forces: FLH FLFt FLFl FRH FRF
% 1 integrated force: int_0^t FLF dt
auxdata.index.variables.states.xy = 1:2;
auxdata.index.variables.states.theta = 3;
auxdata.index.variables.states.uv = 4:5;
auxdata.index.variables.states.omega = 6;
auxdata.index.variables.states.F = 7:11;
auxdata.index.variables.states.intF = 12;

% conventions for sides, given forces or force rates
auxdata.index.side.left.all = 1:3;
auxdata.index.side.left.hind = 1;
auxdata.index.side.left.frontTrail = 2;
auxdata.index.side.left.frontLead = 3;
auxdata.index.side.right.all = 4:5;
auxdata.index.side.right.hind = 4;
auxdata.index.side.right.front = 5;


% Controls:
% 5 Force rate: Fdot
% 2*5 slack variables for work: p, q
% 5 relaxation parameters for limb length constraints
% 1 relaxation parameter for simultaneous limb contact constraint
auxdata.index.variables.control.Fdot = 1:5;
auxdata.index.variables.control.p = 6:10;
auxdata.index.variables.control.q = 11:15;
auxdata.index.variables.control.s = 16:21; % all the relaxation parameters
auxdata.index.variables.control.slimb = 16:20;
auxdata.index.variables.control.ssimult = 21;

% Parameters
% Footfall locations: PLH PTF
auxdata.index.variables.parameter.posLeftHind = 1;
auxdata.index.variables.parameter.posLeftTrailFront = 2;

% Integrals:
% 1 integral for work
% 1 integral for force rate
% 1 integral for slack penalties
auxdata.index.integral.work = 1;
auxdata.index.integral.Fdot = 2;
auxdata.index.integral.slackPen = 3;

% Endpoint constraints
% 5 Kinematic periodicity: yauvw(0.5) - yauvw(0) = 0
% 3 Force continuity: FLH(0.5) - FRH(0) = FRF(0.5) - FTF(0) = FLF(0.5) - FRF(0) = 0
auxdata.index.constraints.endpoint.kinPeriod = 1:5;
auxdata.index.constraints.endpoint.forceContinuity = 6:8;

% Path constraints
% 1 complementary: FTF*int_0^t FLF dt = 0
% 5 complementary: F_i*(l_i - l_imax) >= 0
% 2 normal: y - rH*sin(a) >= 0, y + rF*sin(a) >= 0
% 5 normal: F_i*dl - p_i + q_i = 0
auxdata.index.constraints.path.trailLeadOverlap = 1;
auxdata.index.constraints.path.legLength = 2:6;
auxdata.index.constraints.path.aboveGround = 7:8;
auxdata.index.constraints.path.powerToSlacks = 9:13;

output.result.setup.auxdata = auxdata;
end