% This script uses data collected by CompileOptimalSolutions and determines
% the gait used for each case.
% 
% Data used in the paper are included out of the box as 
% Data/BestGaitSolutions.mat
% 
% If the user wants to compile their own data, run CompileOptimalSolutions
% and then specify the path to the relevant data file on line 15 below

blankSlate

%%%% User Data %%%%
beattol = [0.3 0.03]; % [MinPeakHeight MinPeakDistance]
runtol  = 1e-2; % minimal force for contact (units of mg)
data_path = [pwd,'/Data/BestGaitSolutionsPaper.mat']; % path to data file. 
% ^^ If you run MAIN_SolveSymQuad followed by CompileOptimalSolutions, the
% relevant file will be in Data/mf0.5
show_plot = true; % show a plot of the raw data, and how it looks once missing data are filled
%%%%%%%%%%%%%%%%%%%


% load the data. It can take some time
disp('Data Loading...')
load(data_path);
disp('Data Loaded')

n = length(Optimal);

[beats, isrun] = deal(NaN(n,1));
gait = struct;


%% Determine gait type for each solution
for i = 1:n
    gait(i).out = gaitbeatdetector(Optimal(i).sol.outputBest,beattol);
    beats(i) = gait(i).out.npeaks*2;
    isrun(i) = SymRunDetect(Optimal(i).sol.outputBest,runtol);
end

Ivec = [Optimal(:).I]';
Tvec = [Optimal(:).T]';

gaittype = beats.*isrun;

%% Show a plot of the results
if show_plot
    figure
    scatter(Tvec,Ivec,25,gaittype,'filled')
    box on
    set(gca,'yscale','log')
    colormap('JET')
    colorbar
    xlim([1.4 4.1])
    ylim([0.2 10]+[-0.025 0.8])
    xlabel('T')
    ylabel('I')
end

%% Make vectors of the unique values of T and I

Tval = sort(unique(Tvec));
Ival = sort(unique(Ivec));


%% Create a matrix of T and I, and add relevant data at the correct points

[Tmat, Imat] = meshgrid(Tval,Ival);
[beatmat,isrunmat] = deal(NaN(size(Tmat)));

for i = 1:length(beatmat(:))
    j = Tvec == Tmat(i) & Ivec == Imat(i);
    if any(j)
        beatmat(i) = beats(j);
        isrunmat(i) = isrun(j);
    end
end


% Put the raw data into a matrix
gaittypematraw = beatmat.*isrunmat;

% Fill in missing data
gaittypematfilledraw = inpaint_nans_const(gaittypematraw);

% Switch gait types 2 and -4. This lumps four-beat and two-beat gaits
% together, rather than separating gaits more by walk vs run 
gaittypematnew = gaittypematraw;
i = gaittypematnew(:) == -4;
j = gaittypematnew(:) == 2;
gaittypematnew(i) = 2;
gaittypematnew(j) = -4;

% fill in missing data
gaittypematfillednew = inpaint_nans_const(gaittypematnew);

% show the filled data
if show_plot
    figure
    scatter(Tmat(:),Imat(:),25,gaittypematfilledraw(:))
    colormap('JET')
    set(gca,'yscale','log')
    colorbar
    xlabel('T')
    ylabel('I')
    xlim([1.4 4.1])
    ylim([0.2 10]+[-0.025 0.8])
end

c = strsplit(num2str([beattol,runtol]));
sfx = sprintf('%s_',c{:});

[filepath,~,~] = fileparts(data_path);
savename = ['GaitTypeData',sfx(1:end-1),'.mat'];

save([filepath,filesep,savename],'Tmat','Imat','gaittype*','beattol','runtol','*vec')
