function FoundFieldsList = ...
  SearchFields(Struct, SearchString, SearchOption, StructName)
%
% file SearchFields.m
%
% Created  20061215  Paul A.M. Bune - Alcatel SEL AG
% LastMod  20070326  Paul A.M. Bune - Alcatel-Lucent Deutschland AG
%
% - Searches for field names in a structure array
% - Recursive function
%
% Call options:
%
% SearchFields(AStruct)
% - Displays all fields in structure array AStruct with their
%   complete ("long") names.
%
% SearchFields(AStruct, SearchString)
% - Displays all fields in structure array AStruct of which the name
%   contains SearchString (default = case-insensitive).
%
% SearchFields(AStruct, SearchString, SearchOption)
% - The following options may be indicated by SearchOption:
%   'default'   - Default (case-insensitive search, see above).
%   'case'      - Case-sensitive search.
%   'exact'     - Only those field names are listed which matches the
%                 SearchString exactly (case-insensitive search).
%   'exactcase' - The same as 'exact', but with case-sensitive search.
%   'begin'     - Only those field names are listed of which the beginning
%                 matches the SearchString (case-insensitive search).
%   'begincase' - The same as 'begin', but with case-sensitive search.
%   'end'       - Only those field names are listed of which the end
%                 matches the SearchString (case-insensitive search).
%   'endcase'   - The same as 'end', but with case-sensitive search.
%   If no option is given, 'default' is assumed.
%
% SearchFields(AStruct, SearchString, SearchOption, AStructName)
% - Normally, the name of the structure array AStruct is resolved using
%   the internal MatLab function "inputname". However, this does not work
%   when the name of AStruct contains dots and/or brackets. As a
%   workaround (as well as to enable recursive operation) the name of
%   AStruct can be entered "manually" as a string in AStructName.
%
% FoundFieldsList = SearchFields(AStruct, ...)
% - A list of the found field names is stored as a cellular string array
%   into the output variable FoundFieldsList. If no output variable is
%   indicated, the names are printed as a list on the screen.
%

if nargin > 4
  error('SearchFields: Too many input arguments');
end

if nargin < 1
  error('SearchFields: Too few input arguments');
end

if nargin < 4 % no StructName
  StructName = inputname(1);
  if isempty(StructName)
    disp('SearchFields WARNING:');
    disp('  The name of the structure array could not be resolved.');
    disp('  Note that MatLab routine "inputname" can only resolve names');
    disp('  of variables which do not contain dots and/or brackets.');
    disp('  In the output, the variable name is now left empty');
  end
end

if nargin < 3 % no SearchOption
  SearchOption = 'default';
else
  SearchOption = lower(SearchOption);
end

if nargin < 2 % no SearchString
  SearchString = ''; % Invokes listing of all fields
end

FoundList = {};

if StructName(1) == '>' && ~isstruct(Struct)
  if iscell(Struct)
    for iField = 1:numel(Struct)
      SubscriptVector = ind2sub1(size(Struct), iField);
      CompleteFieldName = [StructName '{' num2str(SubscriptVector(1))];
      for iDim = 2:length(SubscriptVector)
        CompleteFieldName = [CompleteFieldName ',' ...
          num2str(SubscriptVector(iDim))];
      end
      CompleteFieldName = [CompleteFieldName '}'];
      FoundList = [FoundList; ...
        SearchFields(...
        Struct{iField}, ...
        SearchString, ...
        SearchOption, ...
        CompleteFieldName)];
    end
  end
else
  if StructName(1) == '>'
    ListOfFields = fieldnames(Struct);
  else
    ListOfFields = {StructName};
  end
  for iField = 1:length(ListOfFields)
    if StructName(1) == '>'
      CompleteFieldName = [StructName(2:end) '.' ListOfFields{iField}];
    else
      CompleteFieldName = StructName;
    end
    if isempty(SearchString)
      FoundList = [FoundList; CompleteFieldName];
    else
      switch SearchOption
        case 'default'
          if ~isempty(strfind(lower(ListOfFields{iField}), ...
              lower(SearchString)))
            FoundList = [FoundList; CompleteFieldName];
          end
        case 'case'
          if ~isempty(strfind(ListOfFields{iField}, SearchString))
            FoundList = [FoundList; CompleteFieldName];
          end
        case 'exact'
          if strcmpi(ListOfFields{iField}, SearchString)
            FoundList = [FoundList; CompleteFieldName];
          end
        case 'exactcase'
          if strcmp(ListOfFields{iField}, SearchString)
            FoundList = [FoundList; CompleteFieldName];
          end
        case 'begin'
          if strncmpi(ListOfFields{iField}, ...
              SearchString, length(SearchString))
            FoundList = [FoundList; CompleteFieldName];
          end
        case 'begincase'
          if strncmp(ListOfFields{iField}, ...
              SearchString, length(SearchString))
            FoundList = [FoundList; CompleteFieldName];
          end
        case 'end'
          if strncmpi(fliplr(ListOfFields{iField}), ...
              fliplr(SearchString), length(SearchString))
            FoundList = [FoundList; CompleteFieldName];
          end
        case 'endcase'
          if strncmp(fliplr(ListOfFields{iField}), ...
              fliplr(SearchString), length(SearchString))
            FoundList = [FoundList; CompleteFieldName];
          end
        otherwise
          error(...
            ['SearchFields: Unknown search option "' SearchOption '"']);
      end
    end

    if StructName(1) == '>'
      FoundList = [FoundList; ...
        SearchFields(...
        Struct.(ListOfFields{iField}), ...
        SearchString, ...
        SearchOption, ...
        ['>' CompleteFieldName])];
    else
      FoundList = [FoundList; ...
        SearchFields(...
        Struct, ...
        SearchString, ...
        SearchOption, ...
        ['>' CompleteFieldName])];
    end
  end
end

if nargout == 1
  FoundFieldsList = FoundList;
else
  if isempty(FoundList)
    disp('No matching field names found');
  else
    for iField = 1:length(FoundList)
      disp(FoundList{iField})
    end
  end
end

return

function sub = ind2sub1(siz, ind)
%
%   Created 2007-03-26   Paul A.M. Bune
%   LastMod 2007-03-26   Paul A.M. Bune
%
% sub = ind2sub1(siz, ind)
%
%   Multiple subscripts (written into one vector) from single linear index
%   "ind" for an N-dimensional array defined by its array size "siz".
%
%   In principle, the same as MatLab routine "ind2sub", but the result is
%   put into an automatically sized vector rather than into multiple
%   variables.
%
%   Example:
%
%     A    = ones(2, 3, 4);
%     Sub1 = ind2sub1(size(A), 6)
%
%     Sub1 =
%
%          2     3     1
%
%  
if ~isnumeric(ind) || numel(ind) > 1 || ind < 1
  error('"ind" has to be be a single positive integer');
end

if ind > prod(siz)
  error('"ind" exceeds the number of elements defined by "siz"');
end

sub = 1 + mod(ceil(repmat(ind, 1, numel(siz)) ./ ...
  cumprod([1 siz(1:(end-1))])) - 1, siz);

return