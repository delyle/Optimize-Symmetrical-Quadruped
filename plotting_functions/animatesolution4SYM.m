function animatesolution4SYM(GPOPSoutput,savename,showtext,alpha,threshold,bodytype)
% Makes a 4-second animation of the GPOPS solution found by SymQuadOptCtrl
%
% Requires patchline:
% http://www.mathworks.com/matlabcentral/fileexchange/36953-patchline/content/patchline.m
%
% savename: file name for the video. Will be appended with '_animate.mp4'
%
% showtext: display information about the solution.
%   logical vector (length 3): each element corresponds to a different
%       aspect of the information to be displayed outside the plot area
%       ([1 1 1] default) display all elements
%       [1 ~ ~]: Display speed, stride length and forelimb length only
%       [~ 1 ~]: Display nlpinfo, penalty coefficient vector, and
%            cost components (work, force-rate, slack variable penalties)
%       [~ ~ 1]: Display force display threshold
%   'OUTSIDE': Display Speed and MOI text in the outside-middle of the 
%       plotting area
%
% alpha: if true, limb opacity is set by limb force. If false, (default)
%   limbs are shown when limb force exceeds a threshold
%
% threshold: minimum force (mg) at which limbs are shown (if alpha=false)
%
% bodytype: how to display torso (char, case insensitive)
%   'EvenPM', 'Evenpointmass':
%       (default) torso is two point masses at a radius of gyration from
%       the center of mass
%   'PMLIMBS','POINTMASSLIMBS','DISTPM','DISTRIBUTEDPOINTMASS':
%       Torso is displayed as linked point masses at the hips and
%       shoulders, with their relative sizes matching the relative position
%       of the center of mass along the torso
%   'PM','POINTMASS':
%       Torso is displayed as linked point masses at the hips and
%       shoulders, both of the same size.
%

filename = [savename,'_animate'];

narginchk(2,7)

if nargin < 6
    bodytype = 'EvenPM';
    if nargin < 5
        threshold = 0.05;
        if nargin < 4
            alpha = false;
            if nargin < 3
                showtext = logical([1 1 1]);
            end
        end
    end
end

tq = linspace(0,1,120)';

auxdata = GPOPSoutput.result.setup.auxdata;
D = auxdata.D;
[t2,x2,y2,theta2,F2] = SymOutStates2FullCycle(GPOPSoutput);
P = GPOPSoutput.result.interpsolution.parameter;


F3 = interp1(t2,F2,tq,'pchip');
x3 = interp1(t2,x2,tq,'pchip');
y3 = interp1(t2,y2,tq,'pchip');
theta3 = interp1(t2,theta2,tq,'pchip');

% Make new inertial reference frame at average horizontal speed
x0 = tq*D;

% Get new footfall locations
LH = P(1)-x0;
LF1 = P(2)-x0;
LF2 = P(2)+D-x0;
RH1 = P(1)-D/2-x0;
RH2 = P(1)+D/2-x0;
RF1 = P(2)+D/2-x0;
RF2 = P(2)+3*D/2-x0;

mf = auxdata.mf;
m = [mf, 1-mf];
Fmax = max(max(F2));
lmaxF = auxdata.lmax(1);
lmaxH = auxdata.lmax(2);
lmax = max(auxdata.lmax);
F = abs(F3./Fmax);
threshold = threshold/Fmax; % as force is normalized to Fmax, normalize the threshold too
tau = auxdata.tau;
Uh = D*sqrt(tau/lmaxH);
Ihat = auxdata.I*4;
if isfield(auxdata,'c')
    c = auxdata.c;
else
    c = auxdata.c1;
end
I = auxdata.I*4;


r1 = m(2)*[cos(theta3), sin(theta3)];
r2 = -m(1)*[cos(theta3), sin(theta3)];

xs = [r1(:,1)+x3-x0,r2(:,1)+x3-x0]; % x position of shoulders and hips
ys = [r1(:,2)+y3,r2(:,2)+y3]; % y position of shoulders

close all;
figure('Position', [440 378 560 2*lmax*560/3],'color','w')
axes('Position',[.1 .1 .8 .8])
hold on;

switch upper(bodytype)
    case {'EVENPM','EVENPOINTMASS'}
        rgyr = sqrt(I)*[cos(theta3), sin(theta3)]/2; % radius of gyration, relative to half body length.
        xgyr = [-rgyr(:,1)+x3-x0,+rgyr(:,1)+x3-x0];
        ygyr = [-rgyr(:,2)+y3,rgyr(:,2)+y3];
        mrksz = [10 10]; % marker sizes are equal
        mrkst = 's';
    case {'PMLIMBS','POINTMASSLIMBS','DISTPM','DISTRIBUTEDPOINTMASS'}
        xgyr = xs;
        ygyr = ys;
        mrksz = 25*m;
        mrkst = 'o';
    case {'PM','POINTMASS'}
        xgyr = (x3-x0)*[1 1];
        ygyr = [y3,y3];
        mrksz = 18*[1 1];
        mrkst = 'o';
end

yl = [0,1.5*lmax];
xl = [-2.5 2.5];%[-1.5,1.5];

Work = GPOPSoutput.result.solution.phase.integral(1);
Fdot = GPOPSoutput.result.solution.phase.integral(2);
slackPen = GPOPSoutput.result.solution.phase.integral(3);



writerObj = VideoWriter(filename,'MPEG-4');
writerObj.FrameRate = 30;
open(writerObj);

[txt1,txt2,txt3] = deal('');
if strcmpi(showtext,'outside')
    
    txt1 = ['$U''_H = ', num2str(Uh,'%.2f'),',\quad','\hat{I} = ',num2str(Ihat,2),'$'];
    xy_txt = [0.5,1.05;
        0.5,0.9;
        0.9,0.9];
    alignment_txt = {'center','bottom'};
    fontsize_txt = 14;
else
    showtext = logical(showtext);
    xy_txt = [0.1,1.05;
        0.25,1.05;
        0.9,1.05];
    alignment_txt = {'left','bottom'};
    fontsize_txt = 10;
    if showtext(1)
        % This info to be plotted at top-left
        txt1 = ['$U''_H$ = ', num2str(Uh,2),char(10),...
            '$D'' = $',num2str(D),char(10)...
            '$l''_{max} = $', num2str(lmaxF),', ', num2str(lmaxH)];
    end
    if showtext(2)
        % This info to be plotted to the right of the above
        c_text = ['[',regexprep(num2str(c), ' *', ' '),']'];
        txt2 = ['nlpinfo: ',num2str(GPOPSoutput.result.nlpinfo),char(10),...
            'c = ', c_text,...
            char(10),' CoT: Work: ', num2str(Work/D,4),' $\dot{F}^2$: ', num2str(Fdot/D,4), ' slackPen: ' num2str(slackPen,4)];
    end
    if showtext(3)
        txt3 = ['Threshold = $F_{peak}\times$', num2str(threshold),char(10),char(10)];
    end
end
lw = 2;
for i = 1:length(tq)
    cla;
    
    % Parameters
    text(xy_txt(1,1),xy_txt(1,2),txt1,'interpreter','latex','units','normalized','horizontalalignment',alignment_txt{1},'verticalalignment',alignment_txt{2},'fontsize',fontsize_txt);
    text(xy_txt(2,1),xy_txt(2,2),txt2,'interpreter','latex','units','normalized','horizontalalignment',alignment_txt{1},'verticalalignment',alignment_txt{2},'fontsize',fontsize_txt);
    text(xy_txt(3,1),xy_txt(3,2),txt3,...
        'units','normalized','verticalalignment',alignment_txt{2},...
        'horizontalalignment','right','interpreter','latex','fontsize',fontsize_txt);
    % left limbs
    if alpha
        patchline([LF1(i),xs(i,1)],[0,ys(i,1)],'linestyle','--','edgecolor','b','edgealpha',F(i,2),'linewidth',lw);
        patchline([LF2(i),xs(i,1)],[0,ys(i,1)],'linestyle','--','edgecolor','b','edgealpha',F(i,3),'linewidth',lw);
        patchline([LH(i),xs(i,2)],[0,ys(i,2)],'linestyle','--','edgecolor','b','edgealpha',F(i,1),'linewidth',lw);
    else
        if F3(i,1) > threshold
            patchline([LH(i),xs(i,2)],[0,ys(i,2)],'linestyle','--','edgecolor','b','linewidth',lw)
        end
        if F3(i,2) > threshold
            patchline([LF1(i),xs(i,1)],[0,ys(i,1)],'linestyle','--','edgecolor','b','linewidth',lw)
        end
        if F3(i,3) > threshold
            patchline([LF2(i),xs(i,1)],[0,ys(i,1)],'linestyle','--','edgecolor','b','linewidth',lw)
        end
    end
    
    % Torso
    plot(xs(i,:),ys(i,:),'k:','linewidth',lw-0.5)
    plot(xgyr(i,:),ygyr(i,:),'k-','linewidth',lw)
    plot(xgyr(i,1),ygyr(i,1),'ks','linewidth',lw,'markersize',mrksz(1),'marker',mrkst,'markerfacecolor','w')
    plot(xgyr(i,2),ygyr(i,2),'ks','linewidth',lw,'markersize',mrksz(2),'marker',mrkst,'markerfacecolor','w')
    plot(x3(i)-x0(i),y3(i),'kx','markersize',10)
    
    
    % Right Limbs
    if alpha
        patchline([RF1(i),xs(i,1)],[0,ys(i,1)],'linestyle','-','edgecolor','r','edgealpha',F(i,6),'linewidth',lw);
        patchline([RF2(i),xs(i,1)],[0,ys(i,1)],'linestyle','-','edgecolor','r','edgealpha',F(i,7),'linewidth',lw);
        patchline([RH1(i),xs(i,2)],[0,ys(i,2)],'linestyle','-','edgecolor','r','edgealpha',F(i,4),'linewidth',lw);
        patchline([RH2(i),xs(i,2)],[0,ys(i,2)],'linestyle','-','edgecolor','r','edgealpha',F(i,5),'linewidth',lw);
    else
        if F3(i,4) > threshold
            patchline([RH1(i),xs(i,2)],[0,ys(i,2)],'linestyle','-','edgecolor','r','linewidth',lw)
        end
        if F3(i,5) > threshold
            patchline([RH2(i),xs(i,2)],[0,ys(i,2)],'linestyle','-','edgecolor','r','linewidth',lw)
        end
        if F3(i,6) > threshold
            patchline([RF1(i),xs(i,1)],[0,ys(i,1)],'linestyle','-','edgecolor','r','linewidth',lw)
        end
        if F3(i,7) > threshold
            patchline([RF2(i),xs(i,1)],[0,ys(i,1)],'linestyle','-','edgecolor','r','linewidth',lw)
        end
    end
    
    % ground markers
    plot(LH(i),0,'bo')
    
    plot(RH1(i),0,'ro')
    plot(RH2(i),0,'rs')
    
    plot(LF1(i),0,'bo','markerfacecolor','b')
    plot(LF2(i),0,'bs','markerfacecolor','b')
    
    plot(RF1(i),0,'ro','markerfacecolor','r')
    plot(RF2(i),0,'rs','markerfacecolor','r')
    
    
    plot(xl,[0,0],'k:')
    axis equal
    ylim(yl);
    xlim(xl);
    box on
    xlabel('x/l_b')
    ylabel('y/l_b')
    set(gca,'linewidth',1)
    
    drawnow
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
end
close(writerObj)
end