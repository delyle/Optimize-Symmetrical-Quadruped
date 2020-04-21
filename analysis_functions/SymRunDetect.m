function isrun = SymRunDetect(output,mincontactforce)
% Detect whether duty factor exceeds 0.5 in any limb
% In a symmetrical gait over the first half cycle (starting from Left Hind
% touchdown), this can only occur if there is simultaneous contact of the
% Left Hind and Right Hind feet, and/or the Left front and Right front foot
% (in this case we only need to consider left front acting through hind 
% contact). 
% 
% Returns 1 if a run is detected, and -1 otherwise.

X = output.result.interpsolution.phase.state;

F = X(:,7:11); % Axial forces

% 5 forces: FLH FLFt FLFl FRH FRF
FLH = F(:,1); % Axial force of the left hind limb
FLFt = F(:,2);
%FLFl = F(:,3); % don't need this one
FRH = F(:,4);
FRF  = F(:,5);

Hds = FLH > mincontactforce & FRH > mincontactforce; % simultaneous contact in hind
Fds = FLFt> mincontactforce & FRF > mincontactforce; % simultaneous contact in fore


isrun = 1; % it is a run
if any(Hds) || any(Fds)
    isrun = -1; % it is not a run
end