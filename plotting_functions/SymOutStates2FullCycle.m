function [t2,x2,y2,theta2,F2,u2,v2,w2] = SymOutStates2FullCycle(output,soltype)

if nargin < 2
    soltype = 'solution';
end
soltype = lower(soltype);

auxdata = output.result.setup.auxdata;
try
istate = auxdata.index.variables.state;
catch
   istate.F = 7:11;
   istate.xy = 1:2;
   istate.theta = 3;
   istate.dxy = 4:5;
   istate.dtheta = 6;
end

t = output.result.(soltype).phase.time;
X = output.result.(soltype).phase.state;
D = auxdata.D;

% Values are half stride; double up
zs = zeros(length(t),1);
t2 = [t;t(2:end)+0.5];
F = X(:,istate.F);

F1 = [F(:,1:4) zs F(:,5) zs]; % Need to put in this order: [LH LFT LFL RHT RHL RFT RFL];
F2 = [F1;F1(2:end,[4 7 6 5 1 2 3])]; % Add next cycle. Will be repeat but in this order: [RHT RFL RFT RHL LH LFT LFL]

theta = X(:,istate.theta);
theta2 = [theta;theta(2:end)];
x = X(:,istate.xy(1));
y = X(:,istate.xy(2));
x2 = [x;x(2:end)+D/2];
y2 = [y;y(2:end)];

if nargout > 5
    u = X(:,istate.dxy(1));
    v = X(:,istate.dxy(2));
    w = X(:,istate.dtheta);
    u2 = [u;u(2:end)];
    v2 = [v;v(2:end)];
    w2 = [w;w(2:end)];
end
