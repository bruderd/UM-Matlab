function [versionstring, hgid] = version()
% `version()`
%
% Returns the version of MPCTools as a string.
%
% The second output argument returns the HG changeset ID (as a hexidecimal
% string), which is more granular than the versio nnumber.
narginchk(0, 0);
versionstring = '2.1';
hgid = '119b284c77dd';
if nargout() == 0
    fprintf('MPCTools Version %s (hg changeset %s)\n', versionstring, hgid);
end
end%function

