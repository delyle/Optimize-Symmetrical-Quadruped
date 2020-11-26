% This script creates a plot showing gait regions as a function of Murphy
% number and normalized speed.
%
% Out of the box, it should reproduce figure 2a in the paper. Custom data,
% produced by MAIN_DetectGaitTypes, can be specified on line 11 below

blankSlate % clear the worspace

%%%% USER INPUTS %%%% 
FS = 14; % fontsize for axis labels
data_path = 'Data/GaitTypeData0.3_0.03_0.01.mat'; % path to relevant gait data
% ^^ If you have created your own data, by default it will be saved in Data/mf*/
%%%%%%%%%%%%%%%%%%%%%

load(data_path)

Umat = (Tmat/2.4).^(1/-0.32); % this is U/sqrt(g*2L), something like the usual convention

figure('color','w')

TmatFill = [Tmat,Tmat(:,end)]- [diff(Tmat,1,2),(Tmat(:,end)-Tmat(:,end-1))*[1 -1]]/2;
TmatFill = [TmatFill;TmatFill(end,:)];
ImatFill = [Imat;Imat(end,:)] - [(Imat(2,:)-Imat(1,:)); diff(Imat,1,1);-(Imat(end,:)-Imat(end-1,:))]/2;
ImatFill = [ImatFill,ImatFill(:,end)];
UmatFill = (TmatFill/2.4).^(1/-0.32);

% add a row of NaNs to top and right; allows edge cases to be plotted with pcolor or surf
gaittypematFill = [[gaittypematfillednew,zeros(size(gaittypematfillednew,1),1)];zeros(1,size(gaittypematfillednew,2)+1)];
s = pcolor(UmatFill,ImatFill,gaittypematFill);
s.EdgeColor = 'flat';

set(gca,'yscale','log','xscale','log','linewidth',1)
hold on

% Assign the colormap. This was generated with linspecer as 
% cm = linspecer(13,'sequential');
% See https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap
cm = [...
    0.3686    0.3098    0.6353
    0.2005    0.4958    0.7381
    0.3127    0.6971    0.6726
    0.5332    0.8193    0.6451
    0.7616    0.9058    0.6299
    0.9136    0.9604    0.6046
    0.9500    0.9500    0.7116
    0.9938    0.8954    0.5766
    0.9942    0.7583    0.4312
    0.9805    0.5539    0.3044
    0.9288    0.3634    0.2749
    0.8074    0.2075    0.3083
    0.6196    0.0039    0.2588];
colormap(cm)
caxis([-5 5])
xlabel('$U / \sqrt{g2L}$','interpreter','latex','fontsize',FS)
ylabel('$I \,/\, ML^2$','interpreter','latex','fontsize',FS)

% Plot raw data. Make the points slightly darker
hold on
scatter(Umat(:),Imat(:),10,gaittypematnew(:)*1.5,'filled')

% Nice limits to show edges and ticks
xlim([0.20 4.4].*[0.9 1.2])
ylim([0.25 10].*[0.82 1.24])

% Ticks that look good
set(gca,'xtick',[0.2:0.1:1, 2:1:4],'xticklabels',{'0.2','','0.4','','','','','','1','2','','4'})
set(gca,'ytick',[0.3:0.1:1, 2:1:10],'yticklabels',{'','','0.5','','','','','1','2','','','5','','','','','10'})

ax1 = gca;
ax1.Position = ax1.Position - [0 0.01 0 0];
box on

%% Add empirical MOIs. See Table 1
LW = 1;

% dog
I = 0.84;
Urange = xlim; 
plot(Urange,I*[1 1],'k-','linewidth',LW)

% Horse
I = 0.82;
plot(Urange,I*[1 1],'k-','linewidth',LW)

% giraffe
I = 1.66;
plot(Urange,I*[1 1],'k-','linewidth',LW)

% Elephant
I = 0.99;
plot(Urange,I*[1 1],'k-','linewidth',LW)

% Quadrupedal robot (Xi et al. 2016 DOI: 10.1177/0278364915612572)
I = 0.57;
plot(Urange,I*[1 1],'k-','linewidth',LW)



