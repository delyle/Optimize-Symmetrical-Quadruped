function valid = verifySYM(GPOPSoutput,verbose,c)

if nargin <3
    c= false;
end

valid = [true true true true];

[F,~,absF,absL,Ldot] = forceveclengthSYM(GPOPSoutput,'grid');

t = GPOPSoutput.result.solution.phase.time;
U = GPOPSoutput.result.solution.phase.control;
aux = GPOPSoutput.result.setup.auxdata;

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

path2 = p.*q;
crit2 = abs(path2) > 1e-4;
if any(crit2(:))
   valid(2) = false;
   if verbose
       warning('Violation of p.*q == 0')
       disp([find(crit2) path2(crit2)])
   end
end
if c
   % check if energy is conserved
   X = GPOPSoutput.result.solution.phase.state;
   y = X(:,2);
   u = X(:,4);
   v = X(:,5);
   w = X(:,6);

   I = aux.I;
   Fr = aux.Fr;
   E1 = y/Fr + 1/2*(u.^2 + v.^2+ I*w.^2);
   E1 = E1-E1(1); % set initial to 0
   D = aux.D;
   U = aux.Uf;
   T = D/U;
   %E2 = sum(cumtrapz(t,p-q),2);
   [~,E2] = ode45(@(tq,W) workODE(tq,W,t,1/Fr*sum(p-q,2)),t,0);
   disp(max(abs(E1 - E2)));
   close all; plot(t,[E1,E2]);legend('\Delta E','Work')
   
   % now try in SI units
   M = 1;
   G = 10;
   L = 1;
   T = sqrt(1/Fr*L/G);
   
   E1 = M*G*y*L + 1/2*(M*(u*L/T).^2 + M*(v*L/T).^2 + I*M*L^2*(w/T).^2);
   E1 = E1-E1(1); % set initial to 0
   E2 = sum(cumtrapz(t*T,1/Fr*(p-q)*M*(L^2/T^3)),2);
   disp(max(abs(E1-E2)));
   figure; plot(t,[E1,E2]);legend('\Delta E','Work');
   ylabel('joules')
   
   % now try summed forces
   [F,L,~,~,~] = forceveclengthSYM(GPOPSoutput,'grid',true);
   P = dot(F.Htot(:,1:2),L.dot.H,2) + dot(F.Ftot(:,1:2),L.dot.F,2);
   E1 = y/Fr + 1/2*(u.^2 + v.^2+ I*w.^2);
   E1 = E1-E1(1); % set initial to 0
   E2 = cumtrapz(t,1/Fr*P);
   figure; plot(t,[E1,E2]); legend('\Delta E','Work')
   
   median(E1./E2,'omitnan')
   
   % now check dynamics
   
   a_trans = (F.Htot+F.Ftot - ones(length(t),1)*[0 1 0])/Fr;
   mf = aux.mf;
   mh = 1-mf;
   [rf,rh] = deal(zeros(size(F.Htot)));
   theta = X(:,3);
   rf(:,1:2) = mh*[cos(theta),sin(theta)];
   rh(:,1:2) = -mf*[cos(theta),sin(theta)];
   tauF = cross(rf,F.Ftot,2);
   tauH = cross(rh,F.Htot,2);
   
   a_rot = (tauF + tauH)/(I*Fr);
    
   accel = [a_trans+a_rot];
   
   uvw = X(:,4:6);
   uvw_int = cumtrapz(t,accel)+uvw(1,:);
   
   figure;
   subplot(3,1,1)
   plot(t,[uvw(:,1),uvw_int(:,1)]);
   subplot(3,1,2)
   plot(t,[uvw(:,2),uvw_int(:,2)]);
   subplot(3,1,3)
   plot(t,[uvw(:,3),uvw_int(:,3)]);
end
% NetW = sum(trapz(t,abs(W)));
% NetW2 = sum(trapz(t,p + q));


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

