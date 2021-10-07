function [Work, ForceRate, InstData] = CostComponents(output,directWork)
% InstData = instantaneous data

if nargin < 2
   directWork = false; 
end

InstData = struct;

InstData.t = output.result.interpsolution.phase.time;
U = output.result.interpsolution.phase.control;

auxdata = output.result.setup.auxdata;
c = auxdata.c;



InstData.Fdot = c(1)*U(:,auxdata.index.variables.control.Fdot).^2;
p = U(:,auxdata.index.variables.control.p);
q = U(:,auxdata.index.variables.control.q);
if directWork
    [F,~,~,~,Ldot] = forceveclengthSYM(output,'interp');
    P = dot(F,Ldot,2);
    P = permute(P,[1 3 2]);
    InstData.PosWork = sum(c(2)*pospart(P),2);
    InstData.NegWork = sum(c(3)*pospart(-P),2);
else
InstData.PosWork = c(2)*p;
InstData.NegWork = c(3)*q;
end
InstData.pqSlackPen = c(4)*sum(p.*q,2);
InstData.sSlackPen = c(5)*sum(U(:,auxdata.index.variables.control.s),2);

ForceRate = trapz(InstData.t,sum(InstData.Fdot,2));
Work = trapz(InstData.t,sum(InstData.PosWork + InstData.NegWork,2));


end