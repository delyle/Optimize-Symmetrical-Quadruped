function [F,L,absF,absL,Ldot] = forceveclengthSYM(GPOPSoutput,soltype,returnstruct)

if nargin < 3
    returnstruct = false;
    if nargin < 2
        soltype = 'grid';
    end
end

if nargout > 4
    returnderivatives = true;
else
    returnderivatives = false;
end

switch lower(soltype)
    case 'grid'
        X = GPOPSoutput.result.solution.phase.state;
    case 'interp'
        X = GPOPSoutput.result.interpsolution.phase.state;
end
x = X(:,1);
y = X(:,2);
a = X(:,3);

absF = X(:,7:11);
% FLH = F(:,1);
% FTF = F(:,2);
% FLF = F(:,3);
% FRH = F(:,4);
% FRF = F(:,5);

ntime = length(x);
os = ones(ntime,1); zs = zeros(ntime,1);

P = os*GPOPSoutput.result.solution.parameter;

aux = GPOPSoutput.result.setup.auxdata;
mf = aux.mf;
mh = 1-mf;
D = aux.D;

lb_hat =[cos(a), sin(a)];
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

F3 = permute(absF,[1 3 2]);

Fvec(:,1:2,:) = [F3,F3].*L(:,1:2,:)./[absL,absL];

F = Fvec;
absL = permute(absL,[1 3 2]);

if returnderivatives
    uv = X(:,4:5);
    w = X(:,6);
    
    drperp = [-sin(a),cos(a)].*w;
    rHdot = -mf*drperp;
    rFdot = mh*drperp;
    
    Ldot = deal(zeros(ntime,3,5));
    Ldot(:,1:2,[1 4]) = repmat(uv+rHdot,[1 1 2]);
    Ldot(:,1:2,[2 3 5]) = repmat(uv+rFdot,[1 1 3]);
end

if returnstruct
    FHtot = sum(Fvec(:,:,[1 4]),3);
    FFtot = sum(Fvec(:,:,[2 3 5]),3);
    
    F = struct;
    % FLH = F(:,1);
    % FTF = F(:,2);
    % FLF = F(:,3);
    % FRH = F(:,4);
    % FRF = F(:,5);
    
    F.Htot = FHtot;
    F.Ftot = FFtot;
    F.LH = Fvec(:,:,1);
    F.LFt = Fvec(:,:,2);
    F.LFl = Fvec(:,:,3);
    F.RH = Fvec(:,:,4);
    F.RF = Fvec(:,:,5);
    F.abs.LH = absF(:,1);
    F.abs.LFt = absF(:,2);
    F.abs.LFl = absF(:,3);
    F.abs.RH = absF(:,4);
    F.abs.RF = absF(:,5);
    F.abs.Htot = sum(absF(:,[1 4]),2);
    F.abs.Ftot = sum(absF(:,[2 3 5]),2);
    
    Lvec = L;
    L = struct;
    L.LH = Lvec(:,1:2,1);
    L.LFt = Lvec(:,1:2,2);
    L.LFl = Lvec(:,1:2,3);
    L.RH = Lvec(:,1:2,4);
    L.RF = Lvec(:,1:2,5);
    L.abs.LH = absL(:,1);
    L.abs.LFt = absL(:,2);
    L.abs.LFl = absL(:,3);
    L.abs.RH = absL(:,4);
    L.abs.RF = absL(:,5);
    if returnderivatives
        L.dot.H = Ldot(:,1:2,1);
        L.dot.F = Ldot(:,1:2,2);
    end
end



