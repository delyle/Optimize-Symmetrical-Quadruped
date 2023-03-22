function GPOPSoutput = updateFineMesh(GPOPSinput,PM_NW,plotfigures)
% This function takes GPOPSinput, resulting from SymQuadOptCtrl, and
% interpolates additional collocation points midway between existing
% points. This serves as a guess for a new round of optimization using
% GPOPS-II, resulting in GPOPSoutput.
%
% PM_NW is a logical vector of length two. First element is point mass,
% second is whether optimization is on net work. If length == 1, we will
% assume we're not optimizing on net work

if nargin < 3
    plotfigures = false;
end


pointmass = PM_NW(1);
if length(PM_NW) == 1
    NetWork = false;
else
    NetWork = PM_NW(2);
end

F = GPOPSinput.result.solution.phase.state(:,7:11);
Fi = GPOPSinput.result.interpsolution.phase.state(:,7:11);
t = GPOPSinput.result.solution.phase.time;
ti = GPOPSinput.result.interpsolution.phase.time;
t_bw = t(1:end-1) + diff(t)/2; % intermediate points
t_mat = [[NaN;t_bw],t]'; % matrix with intermediate points, regular points
tq = t_mat(2:end)';

if plotfigures
    figure;
    plot(t,F,'o');
    set(gca,'colororderindex',1);
    hold on
    plot(ti,Fi);
    
    Fq = interp1(ti,Fi,tq,'linear');
    set(gca,'colororderindex',1);
    plot(tq,Fq,'.')
    drawnow
end

intfun = @(xi) interp1(ti,xi,tq);

Xi = GPOPSinput.result.interpsolution.phase.state;
Ui = GPOPSinput.result.interpsolution.phase.control;

if pointmass
   Xi(:,[3,6]) = []; % remove angular positions if it's a pointmass
   if NetWork
       GPOPSfun = @SymQuadOptCtrl_pointmass_NCW;
   else
        GPOPSfun = @SymQuadOptCtrl_pointmass;
   end
else
    if NetWork
        GPOPSfun = @SymQuadOptCtrl_NCW;
    else
        GPOPSfun = @SymQuadOptCtrl;
    end
end

old_colpoints = GPOPSinput.result.setup.mesh.phase.colpoints;
new_colpoints = old_colpoints*2;

input = GPOPSinput;
input.result.solution.phase.time = tq;
input.result.solution.phase.state = intfun(Xi);
input.result.solution.phase.control = intfun(Ui);
input.result.setup.mesh.phase.colpoints = new_colpoints;

auxdata = GPOPSinput.result.setup.auxdata;
auxdata.meshMaxIter = 8;
if isfield(auxdata,'LimbWork')
    auxdata.LimbWork = logical(auxdata.LimbWork);
end
auxdata.downsample = false;

GPOPSoutput = GPOPSfun(auxdata,input);