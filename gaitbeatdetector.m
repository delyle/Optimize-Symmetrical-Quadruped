function beats = gaitbeatdetector(GPOPSoutput,minpeakdiff,showplot)

if nargin < 3
    showplot = false;
    if nargin < 2
        minpeakdiff = [0.2,0.1];
    end
end

t = GPOPSoutput.result.interpsolution.phase.time; % get time
U = GPOPSoutput.result.interpsolution.phase.control;
p = U(:,6:10); % Positive part of actuator power
q = U(:,11:15); % Negative part of actuator power

minpeakheight = minpeakdiff(1);
minpeakdistance = minpeakdiff(2);

Ptot = sum(p-q,2);
Ptotn = -Ptot/min(Ptot(1:end-1));
Ptotn_short = Ptotn(1:end-1);
Ptotn_ext = repmat(Ptotn_short,[3 1]);
t_ext = [t(1:end-1)-0.5;t;t(2:end-1)+0.5];
[~,i_max] = findpeaks(-Ptotn_ext,t_ext,'minpeakheight',minpeakheight,'minpeakdistance',minpeakdistance);
if showplot
    figure;
    findpeaks(-Ptotn_ext,t_ext,'minpeakheight',minpeakheight,'minpeakdistance',minpeakdistance);
end

istart = length(Ptot);
iend   = 2*(length(Ptot)-1); % we don't want to count t=0 and 0.5; they are both the same point. We search up to 0.5-i

i_max = find(ismember(t_ext,i_max)); % convert from time domain to indices

ii_keep = i_max >= istart & i_max <= iend;
i_max_keep = i_max(ii_keep)-istart+1;% keep only indices within original signal

beats.npeaks = length(i_max_keep);
beats.PeakP = Ptotn(i_max_keep);
beats.ipeaks = i_max_keep;
beats.tpeaks = t_ext(i_max_keep);

end