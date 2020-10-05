function GPOPSoutput = updateFineMesh(GPOPSinput,plotfigures)
% This function takes GPOPSinput, resulting from SymQuadOptCtrl, and
% interpolates additional collocation points midway between existing
% points. This serves as a guess for a new round of optimization using
% GPOPS-II, resulting in GPOPSoutput.

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
    resetcolor
    hold on
    plot(ti,Fi);
    
    Fq = interp1(ti,Fi,tq,'linear');
    resetcolor
    plot(tq,Fq,'.')
    drawnow
end

intfun = @(xi) interp1(ti,xi,tq);

Xi = GPOPSinput.result.interpsolution.phase.state;
Ui = GPOPSinput.result.interpsolution.phase.control;



old_colpoints = GPOPSinput.result.setup.mesh.phase.colpoints;
new_colpoints = old_colpoints*2;

input = GPOPSinput;
input.result.solution.phase.time = tq;
input.result.solution.phase.state = intfun(Xi);
input.result.solution.phase.control = intfun(Ui);
input.result.setup.mesh.phase.colpoints = new_colpoints;

auxdata = GPOPSinput.result.setup.auxdata;
auxdata.meshMaxIter = 8;
auxdata.LimbWork = logical(auxdata.LimbWork);
auxdata.downsample = false;

GPOPSoutput = SymQuadOptCtrl(auxdata,input);

