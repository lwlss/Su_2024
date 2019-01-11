function struct2var(s)

%STRUCT2VAR Convert structure array to workspace variables.
%   STRUCT2VAR(S) converts the M-by-N structure S (with P fields)
%   into P variables defined by fieldnames with dimensions M-by-N.  P
%   variables are placed in the calling workspace.

if nargin < 1
    error('struct2var:invalid','No input structure')
elseif nargin > 1
    error('struct2var:invalidt','Too many inputs')
elseif ~isstruct(s)
    error('struct2var:invalid','Input needs to be a structure data type')
end

[r,c] = size(s);
names = fieldnames(s);

for i=1:length(names)
    assignin('caller',names{i},s.(names{i}))
end