function All_results_plotSYM(GPOPSoutput)
% Plot results from a single solution of SymQuadOptCtrl

t = GPOPSoutput.result.solution.phase.time;
X = GPOPSoutput.result.solution.phase.state;
P = GPOPSoutput.result.solution.parameter;

auxdata = GPOPSoutput.result.setup.auxdata;
F = X(:,7:11);

Uf = auxdata.Uf;
lmax = auxdata.lmax(1);
mf = auxdata.mf;
D = auxdata.D;
Fmax = auxdata.Fmax;



if isfield(auxdata,'c')
    c = auxdata.c;
else
    c = auxdata.c1;
end

figure('position',[782     1   657   804],'color','w');
c_text = ['[',regexprep(num2str(c), ' *', ' '),']'];
str = ['$Fr_A\equiv D/\sqrt{gl_{Fmax}T^2}=$ ',num2str(Uf), ', $l_{Fmax}=',num2str(lmax),...
        'l_b$, $m_F=',num2str(mf),'$, $D=',num2str(D),'l_b$, $F_{max}=',num2str(Fmax),'mg$, $c=$',c_text];

subplot(5,1,1)
    [hAx,~,~] = plotyy(t,X(:,1)-auxdata.D*t,t,X(:,2));
    ylabel(hAx(1),'$X/l_b - Dt/T$','interpreter','latex')
    ylabel(hAx(2),'$Y/l_b$','interpreter','latex')
    text(0.5,1.17,str,'units','normalized','horizontalalignment','center',...
        'interpreter','latex');
subplot(5,1,2)
    plot(t,[X(:,4)-D,X(:,5)])
    legend({'$\dot{X}T/l_b - D/T$','$\dot{Y}T/l_b$'},'interpreter','latex')  
subplot(5,1,3)
    [hAx,~,~] = plotyy(t,X(:,3)*180/pi,t,X(:,6)*180/pi);
    ylabel(hAx(1),'$\theta$ [$^\circ$]','interpreter','latex')
    ylabel(hAx(2),'$\dot{\theta}T$ [$^\circ$]','interpreter','latex')
hs4 = subplot(5,1,4);
    h = plot(t,F(:,[1,2,4,5]),'linewidth',1);
    hold on
    h(5,:) = plot(0,0,'k-','linewidth',1);
    h(6,:) = plot(0,0,'k--','linewidth',1);
    legend(h,{'LH','LF','RH','RF','Trail','Lead'},'location','northoutside','orientation','horizontal')
    ax = gca;
    ax.ColorOrderIndex = 2;
    plot(t,F(:,3),'--')
    ylabel('$F/(mg)$','interpreter','latex')
    xlabel('t^*')
hs5 = subplot(5,1,5);
       plot(P(1),0,'s');
       hold on
       plot(P(2),0,'s');
       plot(P(1)-D/2,0,'s')
       plot(P(2)+D/2,0,'s');
       c = get(gca,'colororder');
       plot(P(2)+D,0,'o','color',c(2,:))
       
    xl = xlim;
    plot(xl,[0 0],'k-')
    h = plot(-2,-2,'ks');
    h.MarkerFaceColor = h.Color;
    h(2,:) = plot(-2,-2,'ks');
    legend(h,{'Lead','Trail'})
    
    % make x ticks and tick labels in terms of D
    xt = ceil(xl(1)/D)*D:D:floor(xl(2)/D)*D;
    xtl = textscan(num2str(xt/D),'%s')';
    xtl = strcat(xtl{1},'D');
    xtl = strrep(xtl,'0D','0');
    set(hs5,'xtick',xt,'xticklabel',xtl,'ytick',[])
    axis equal
    xlim(xl);
    
    title('Footfall positions')
    
hs5.Position = [0.1239    0.0541    0.7750    0.1];
hs4.Position = [0.1300    0.2239    0.7750    0.1467];