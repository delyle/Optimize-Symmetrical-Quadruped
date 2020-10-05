function valid = verifySYM(GPOPSoutput,verbose)
% Checks if complementarity conditions are satisfied
%
% Complementarity conditions are imposed as penalties to the (augmented) 
% cost function. Therefore, it is important to ensure that they are 
% satisfied post-hoc

% Initialize validation vector.
valid = [true true true true];

% Get relevant data from the simulation

% get force and leg length as vectors
[F,~,absF,absL,Ldot] = forceveclengthSYM(GPOPSoutput,'grid');

% Get time, controls, auxdata
t = GPOPSoutput.result.solution.phase.time;
U = GPOPSoutput.result.solution.phase.control;
aux = GPOPSoutput.result.setup.auxdata;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that F*ldot = p - q (the slack variables p and q are the positive 
% and negative parts of limb work)

p = U(:,6:10);
q = U(:,11:15);

if ~aux.LimbWork
    uv = GPOPSoutput.result.solution.phase.state(:,4:5);
    Ldot = zeros(size(F));
    Ldot(:,1:2,:) = uv.*ones(length(t),2,5);
end
P = dot(F,Ldot,2);
P = permute(P,[1 3 2]);
path1 = abs(P - p + q);

crit1 = path1(1:end-1,:) > 1e-4;
if any(any(crit1))
   valid(1) = false;
   if verbose
   warning('Violation of P - p + q')
   disp([find(crit1) path1(crit1)])
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that p*q = 0 is satisfied

path2 = p.*q;
crit2 = abs(path2) > 1e-4;
if any(crit2(:))
   valid(2) = false;
   if verbose
       warning('Violation of p.*q == 0')
       disp([find(crit2) path2(crit2)])
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that limbs do not produce force when the limb length exceeds
% maximum value, i.e. F*(lmax - l) >= 0

forcetol = 0.02;
limboff = absF < forcetol;

absLon = absL;
absLon(limboff) = NaN;

aux = GPOPSoutput.result.setup.auxdata;

leg_tol = 0.01;
leg_mult = 1+leg_tol;
crit3 = [absLon(:,[1 4]) > aux.lmax(2)*leg_mult, absLon(:,[2 3 5]) > aux.lmax(1)*leg_mult];
if any(any(crit3))
    valid(3) = false;
    if verbose
    warning('Violation of F*(lmax - l) >= 0')
    disp([find(crit3) absL(crit3)])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check that limb acting through the trailing footfall is not active when
% limb acting through the leading footfall is active (For left-front limb)

F_LFt = absF(:,2);
intF_LFl = GPOPSoutput.result.solution.phase.state(:,12);

path3 = abs(F_LFt.*intF_LFl);
crit4 = path3 > 0.01^2;
if any(crit4)
   valid(4) = false;
   if verbose
   warning('Violation of F_LFt.*intF_LFl')
   disp([find(crit4) path3(crit4)])
   end
end

end

