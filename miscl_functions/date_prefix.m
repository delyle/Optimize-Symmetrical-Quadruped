function str = date_prefix(timecode)
% DATEPREFIX this function returns a string containing the current date in
% a specified format
% timecode = time specifier, see > help datestr
%          = 'yyyymmdd' (default)
narginchk(0,1)
if nargin == 0
    timecode = 'yyyymmdd';
end

str = datestr(datetime('now'),timecode);

end