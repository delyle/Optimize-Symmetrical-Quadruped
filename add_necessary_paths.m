%add paths

addpath(pwd)
addpath([pwd,'/plotting_functions'])
addpath([pwd,'/analysis_functions'])

% Temporary: add GPOS-II, SNOPT paths

% Add SNOPT
addpath([userpath, '/SNOPT/snopt7-6/matlab'])

% Add SNOPT license
setenv('SNOPT_LICENSE','/Users/delyle/Documents/MATLAB/SNOPT/license/snopt7.lic')

% Add GPOPS-II paths
gpopsdir = [userpath,'/GPOPS-II'];
gpopssubdir = {'gpopsUtilities','lib','license','nlp'};
addpath(gpopsdir)
for i = 1:length(gpopssubdir)
    addpath(genpath([gpopsdir,'/',gpopssubdir{i}]))
end