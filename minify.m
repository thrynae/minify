function str=minify(str,BeastMode,max_length_after_merge,EntryFunctionNames,end_function_char)
%Process Matlab code into a solid block of compact and unreadable but functionally equivalent code.
%The output may require some manual tweaking to be optimally compact.
%
%Syntax:
% str_out=minify(str_in)
% str_out=minify(str_in,BeastMode)
% str_out=minify(str_in,BeastMode,max_length_after_merge)
% str_out=minify(str_in,BeastMode,max_length_after_merge,EntryFunctionNames)
%
%Arguments:
% str_out:                mx1 cell array with the minified code
% str_in:                 nx1 cell array of char arrays
% BeastMode:              boolean that triggers additional steps to further minify the code
%                         (specific aspects can be set individually by supplying a struct instead)
% max_length_after_merge: target maximum length (might be exceeded by original code after the
%                         removal of the line continuation)
%                         if <10 this is multiplied by the maximum line length preference
%                         (if you want more than 10x, you can provide a negative value)
% EntryFunctionNames:     cell array with function names that should be kept unchanged, these
%                         functions will be move to the top if there are multiple functions. This
%                         defaults to the first function (if any).
%
%BeastMode parameters:
% (missing fields will be filled with the default, BeastMode=false will set all fields to false)
% trim_spaces                   [true]  Strip leading spaces and remove all double spaces that do
%                                       not occur inside a char or a string.
% compress_to_block             [false] Attempt to compress the entire input to a single block.
%                                       Setting this to true will cause the resulting code to be
%                                       incompatible between modern Matlab releases and Octave and
%                                       ML6.5. The reason is that the former require a space
%                                       between 'end' and 'function', while the latter require a
%                                       comma or semicolon.
% compress_functions_separately [true]  Compress each function to a separate block of code.
%                                       If proper detection of nested function is implemented, this
%                                       will compress each parent function including the nested
%                                       functions to a single block.
%                                       If compress_to_block is set to true, this is ignored.
% keep_original_function_names  [false] Do not rename local functions.
%                                       This setting can be used if any code uses eval with the
%                                       name of a local function in a char or a string. (you can
%                                       use func2str(@local_function) to work around this issue)
% contains_nested_functions     [false] Setting this to true will skip the steps that would break
%                                       nested functions.
%
%Example:
%  str=readfile('readfile.m');
%  str=minify(str);
%  fid=fopen(fn_out,'w');fprintf(fid,'%s\n',str{:});fclose(fid);
%
% The likely use cases fall generally under two categories:
%    1) Attaching code as a dependency to allow your code to run, without having to refer to
%       separate FEX/Github entries. In some such places understanding the attached code is not
%       important, while space is at a premium (e.g. pdf attachments, where large amounts of
%       dependent code would be distracting).
%    2) Code obfuscation. Since p-code will only run on a subset of Matlab releases (and not at all
%       on GNU Octave), using that will limit compatibility. Additionally, given that the
%       encryption has been broken, this might be a good additional step (or replacement) to hide
%       the function of your code without harming its function.
%
% The process of compacting Matlab code is split into several steps:
%    1) Strip all comments and line continuations.
%    2) Sort the functions (if there are multiple) in order to keep the entry function(s) at the
%       top of the output.
%    3) Separate the code and the embedded chars/strings.
%    4) Parse the code to select only the places where a new variable may be created in the
%       workspace. False positives should be avoided at all costs. False negatives should be
%       avoided, but are not a major issue.
%    5) Use that (and the list of local functions) to create a dictionary of variable and function
%       names.
%    6) Remove the entry function from that dictionary and replace every occurrence of a variable
%       or local function with a shorter/pseudonymized one.
%    7) (BeastMode==true) Replace all double spaces in code with single space.
%    8) Put the code and chars/strings back together.
%    9) (BeastMode==true) Merge lines if the result is shorter that a fixed length. A space, comma,
%       or semicolon may be added. Be careful with functions that rely on printing results to the
%       command window by omitting the semicolon, although such functions should probably not be
%       minified in the first place (as the variable names are changed).
% Steps 3-8 should be performed for each function separately.
%
% This function was tested on a random sample of 1000 m files from the FileExchange. Some limits
% were imposed on the selection of files: only submissions with 5 downloads or more in the last 30
% days, at most 5 files per submission (taking the first files in the hierarchy, without checking
% if those actually would contain the main function), ignoring functions with fewer than 5 or more
% than 1000 lines, and ignoring functions with lines over 1000 characters long (as those probably
% contain data, not real code). If BeastMode is set to true the size of functions tends to be
% reduced to about 13% of the original number of lines (half of them are between 9.5% and 16.7%).
% Setting BeastMode to false will increase the variability by a lot, and results in file sizes of
% about 51% of the original number of lines (38.3-66.5%).
% The compression depends mostly on the amount of comments, typical line length, amount and length
% of chars/strings, and typical variable name length.
%
%Compatibility considerations:
% - Support for eval and friends is limited to situation where they don't rely on variable/function
%   names (e.g. if it is used to create an anonymous function). You could use something like
%   func2str(@local_function) to use a local function call inside an eval statement.
% - Nested functions are tricky to extract. It is possible (match up end statements with 'if',
%   'try', 'while', 'for', 'parfor', and 'function', then confirm every function has an 'end' and
%   ignore nested functions while sorting the functions). I don't use them (as they are
%   incompatible with Matlab 6.5 and Octave), but it isn't impossible to modify this function. You
%   can set contains_nested_functions to true to turn off the parts that will interfere with nested
%   functions.
% - There is no support for the arguments block. This is not a fundamental issue, it is just not
%   yet implemented in this function.
%
%  _____________________________________________________________________________
% | Compatibility   | Windows XP/7/10 | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |-----------------|-----------------|------------------|----------------------|
% | ML R2020b       | W10: works      |  not tested      |  not tested          |
% | ML R2018a       | W10: works      |  works           |  not tested          |
% | ML R2015a       | W10: works      |  works           |  not tested          |
% | ML R2011a       | W10: works      |  works           |  not tested          |
% | ML R2010b       | not tested      |  works           |  not tested          |
% | ML R2010a       | W7:  works      |  not tested      |  not tested          |
% | ML 7.1 (R14SP3) | XP:  works      |  not tested      |  not tested          |
% | ML 6.5 (R13)    | W10: works      |  not tested      |  not tested          |
% | Octave 6.1.0    | W10: works      |  not tested      |  not tested          |
% | Octave 5.2.0    | W10: works      |  works           |  not tested          |
% | Octave 4.4.1    | W10: works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
% note: Octave and ML6.5 require a semicolon or comma between 'end' and 'function', while newer
% Matlab releases either require or allow a space. Setting compress_to_block to true will make the
% resulting code incompatible between the two styles.
%
% Version: 1.0
% Date:    2020-12-23
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email = 'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

% Set defaults.
persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
if nargin<1,error('not enough inputs');end
if nargin<2 || isempty(BeastMode),BeastMode=true;end
if nargin<3 || isempty(max_length_after_merge),max_length_after_merge=1.5;end
%(strip the comments before setting the default for input arg 4)
if nargin<5,end_function_char='';end % Undocumented, required if compress_to_block==true.

% Check inputs and parse the options.
valid=false;try valid=isequal(str,cellstr(str));catch,end
if ~valid,error('first input must be a cellstr');end
BeastMode=FillBeastModeStruct(BeastMode);
if max_length_after_merge<10 % Multiply with the maximum line length preference.
    %To have more than 10x the line length, you can provide a negative value.
    max_length_after_merge=ceil(abs(max_length_after_merge)*GetMaxLineLengthPref);
end
% This is undocumented on purpose. If you set compress_to_block to true, this can be used this
% input to 'cross-compile' between ML6.5 and Octave, and modern releases of Matlab. I.e. you can
% use ML6.5 to create code that will run only on newer releases or use a newer release to create
% code that will only run on ML6.5 and Octave.
% Although this is meant to be either a space or a comma, you could also use char(10) to get
% approximately the same result as with compress_functions_separately=true.
SetDefault_end_function_char(end_function_char);

% Strip all comments and merge line continuations (i.e. the ...).
str=StripComments(str);

% If there are only comments, return an empty array.
if numel(str)==0,return,end

% Only now can we reliably detect the name of the first function.
if nargin<4 || isempty(EntryFunctionNames),EntryFunctionNames=ListLocalFunctions(str{1});end

EntryFunctionNames=cellstr(EntryFunctionNames);% Ensure EntryFunctionNames is a cellstr.

% Sort entry functions to the top, unless there are nested functions.
if ~BeastMode.contains_nested_functions
    functions=SortFunctions(str,EntryFunctionNames);
else
    % As long as the code in this else section makes sure not to separate nested functions from
    % their parent function, the rest of the BeastMode settings should work normally.
    functions={str};
end

% Create the replacement dictionary for the local functions.
if BeastMode.keep_original_function_names
    names__fun=cell(0,2);
else
    [ignore,only_cmd]=SplitCodeAndChar(vertcat(functions{:})); %#ok<ASGLU>
    names__fun=ListLocalFunctions(only_cmd);
    names__fun(ismember(names__fun,EntryFunctionNames))=[];
    for n=1:numel(names__fun),names__fun{n,2}=sprintf('f%02d',n-1);end
    
    if BeastMode.contains_nested_functions
        % Create replacement names.
        names__var___all=ListVariables(only_cmd);
        for n=1:numel(names__var___all),names__var___all{n,2}=sprintf('v%03d',n-1);end
    end
end

%The look forward and look backward doesn't work on ML6.5.
persistent UseLookaroundRegExp
if isempty(UseLookaroundRegExp),UseLookaroundRegExp=ifversion('>=',7,'Octave','>',0);end

for fun_ind=1:numel(functions)
    str=functions{fun_ind};
    
    % Separate code and chars and determine the names of internal functions and variables.
    [full_file,only_cmd]=SplitCodeAndChar(str);
    names__var=ListVariables(only_cmd);
    
    if BeastMode.contains_nested_functions
        % Create a dictionary with only the variable names that occur in this function. This will
        % keep variables in nested functions working as intended, while avoiding renaming function
        % calls with variables (e.g. 'line' might be used in one function as a variable, but as the
        % the graphics primitive in another).
        % This will reduce the code obfuscation and increase the likelyhood of overflowing the 3
        % digit size of variable names, which is why it not the default.
        L=ismember(names__var___all(:,1),names__var);
        dict=[names__fun;names__var___all(L,:)];
    else
        % Create a new dictionary for each function. This means every function will start counting
        % variables at v000.
        for n=1:numel(names__var),names__var{n,2}=sprintf('v%03d',n-1);end
        dict=[names__fun;names__var];
    end
    if numel(dict)==0,dict=cell(0);end% Set size(dict,1) to 0.
    
    % Pad scalar input so find works as expected.
    scalar_input=size(full_file,1)==1;
    if scalar_input,full_file(2,1)={' '};end
    
    % Perform the replacement element by element.
    X=~cellfun('isempty',full_file);
    X(:,2:2:end)=false;% Don't replace things inside of strings/chars.
    if UseLookaroundRegExp
        for d=1:size(dict,1)
            RE=['(^|(?<=[^a-zA-Z0-9_\.]))(' dict{d,1} ')(?=[^a-zA-Z0-9_]|$)'];
            for x=find(X).'
                full_file{x}=regexprep(full_file{x},RE,dict{d,2});
            end
        end
    else
        for d=1:size(dict,1)
            RE=['([^a-zA-Z0-9_\.])(' dict{d,1} ')([^a-zA-Z0-9_])'];
            for x=find(X).'
                str=[' ' full_file{x} ' '];
                % Perform the replacement twice in case of overlapping tokens (e.g. 'a=a+1;').
                str=regexprep(str,RE,['$1' dict{d,2} '$3'],'tokenize');
                str=regexprep(str,RE,['$1' dict{d,2} '$3'],'tokenize');
                full_file{x}=str(2:(end-1));
            end
        end
    end
    
    if BeastMode.trim_spaces
        % Replace double spaces with a single space, but not in the chars.
        for x=find(X).'
            removed_double_spaces=inf;
            while removed_double_spaces~=0
                length1=length(full_file{x});
                full_file{x}=strrep(full_file{x},'  ',' ');
                length2=length(full_file{x});
                removed_double_spaces=length1-length2;
            end
        end
        % Strip leading space as well.
        for r=1:size(X,1)
            if numel(full_file{r,1})>=1 && strcmp(full_file{r,1}(1),' '),full_file{r,1}(1)='';end
        end
    end
    
    % Un-pad scalar input.
    if scalar_input,full_file(2,:)=[];end
    
    % Reshape back to cellstr.
    for r=1:size(full_file,1)
        line=full_file(r,:);
        if isOctave,line(cellfun('isempty',line))={''};end
        full_file{r,1}=[line{:}];
    end
    str=full_file(:,1);
    functions{fun_ind}=str;
end

% Note: Manual tweaking may still be needed to get to the most compact block possible.
if BeastMode.compress_to_block
    % Merge all lines into a single line (adding a comma, semicolon or space when needed) and
    % attempt to split that along commas and semicolons (adding ellipses where needed).
    if numel(EntryFunctionNames)>1
        % Split the entry functions from the rest and compress each to a block.
        funs=cell(1,1+numel(EntryFunctionNames));
        funs(1:numel(EntryFunctionNames))=functions(1:numel(EntryFunctionNames));
        if numel(functions)>numel(EntryFunctionNames)
            remainder=functions((numel(EntryFunctionNames)+1):end);
            funs{end}=vertcat(remainder{:}); % Add the remainder of the code.
        else
            funs(end)=[]; % Remove the last cell if it is not needed.
        end
        for n=1:numel(funs)
            funs{n}=CompressFunctionBody(funs{n},max_length_after_merge);
        end
        % Merge back to a single array.
        str=vertcat(funs{:});
    else
        % Merge back to a single array.
        str=vertcat(functions{:});
        str=CompressFunctionBody(str,max_length_after_merge);
    end
elseif BeastMode.compress_functions_separately
    % Compress each function separately.
    for n=1:numel(functions)
        if numel(functions{n})==1,continue,end % Skip the compression if it is already 1 line.
        functions{n}=CompressFunctionBody(functions{n},max_length_after_merge);
    end
    % Merge back to a single array.
    str=vertcat(functions{:});
else
    % No BeastMode option here, so merge back to a single array.
    str=vertcat(functions{:});
end

if nargout==0
    clear str
end
end
function opts=FillBeastModeStruct(BeastMode)
% Parse the input to a struct.
% If the input is true, the default is selected, otherwise all fields are set to false.
% Missing fields are filled from the default.
persistent default all_false
if isempty(default)
    default=struct;
    default.trim_spaces=true;
    default.compress_to_block=false;
    default.compress_functions_separately=true;
    default.keep_original_function_names=false;
    default.contains_nested_functions=false;
    
    all_false=default;
    fn=fieldnames(all_false);
    for n=1:numel(fn)
        all_false.(fn{n})=false;
    end
end
if ~isa(BeastMode,'struct')
    [isScalarLogical,BeastMode]=test_if_scalar_logical(BeastMode);
    if ~isScalarLogical || BeastMode
        opts=default;return
    else
        opts=all_false;return
    end
end
opts=default;
fn=fieldnames(BeastMode);
for n=1:numel(fn)
    % Parse all struct fields, but make sure all fields are scalar logicals.
    % Any invalid field content will be converted to true.
    [isScalarLogical,val]=test_if_scalar_logical(BeastMode.(fn{n}));
    if isScalarLogical,opts.(fn{n})=val;else,opts.(fn{n})=true;end
end
end
function out=CompressFunctionBody(body,max_length_after_merge)
%Compress the function body. Semicolons, commas, spaces, and ellipses will be added if needed.

body=MergeShortLines(body,inf,false);
txt=body{1};
body=SplitLineToCodeAndChar(txt);
% Find the places where we can put a line break with ellipsis (after a comma or space) or without
% one (after a semicolon).
loc=cell(size(body));
for m=1:size(loc,2)
    c=zeros(size(body{1,m}));
    c(body{1,m}==';')=1;
    c(body{1,m}==',')=2;
    c(body{1,m}==' ')=3;
    c(body{1,m}=='+')=4;
    c(body{1,m}=='-')=4;
    c(body{1,m}=='@')=4;
    loc{1,m}=c;
    loc{2,m}=zeros(size(body{2,m}));
end
loc=horzcat(loc{:});
flipped_loc=fliplr(loc);
flipped_txt=fliplr(txt);
out=cell(ceil(numel(loc)/max_length_after_merge),1);
n_out=1;
while ~isempty(flipped_txt)
    if n_out>1 && strcmp(out{end-1},'...')
        % Prevent the same split to be performed over and over again. This will occur if a line
        % ends with ',...' or ' ),...' with the semicolon being more than max_length_after_merge
        % from the end.
        % This code will result in a line that is longer than max_length_after_merge.
        n_out=n_out-1;
        idx=find(flipped_loc>0);
        flipped_loc(idx(1))=0;
    end
    % First try to find a semicolon or comma.
    idx=find(flipped_loc>0);
    if isempty(idx) || numel(flipped_txt)<=max_length_after_merge
        % Either this line was short enough in the original file, or we reached the first line.
        out{n_out}=flipped_txt(end:-1:1);
        break
    end
    % Find the last split location before the line limit.
    idx(idx>max_length_after_merge)=[];
    
    if isempty(idx)
        idx=find(flipped_loc>0);idx=min(idx); %#ok<MXFND>
    else
        % Prefer an earlier semicolon if it is fewer than 5 characters before the last comma. For
        % spaces use a threshold of 8. For the other characters use a threshold of 15.
        last_separator(1)=max([-inf idx(flipped_loc(idx)==1)]);
        last_separator(2)=max([-inf idx(flipped_loc(idx)==2)]);
        last_separator(3)=max([-inf idx(flipped_loc(idx)==3)]);
        last_separator(4)=max([-inf idx(flipped_loc(idx)==4)]);
        [ignore,type]=max(last_separator-[0 5 8 15]); %#ok<ASGLU>
        idx=last_separator(type);
        add_ellipsis=flipped_loc(idx)~=1;
    end
    idx=max(1,idx-1);% Don't include the comma or semicolon itself.
    out{n_out}=fliplr(flipped_txt(1:idx));
    n_out=n_out+1;
    flipped_txt(1:idx)='';flipped_loc(1:idx)='';
    if add_ellipsis
        flipped_txt=['...' flipped_txt]; %#ok<AGROW>
        flipped_loc=[0 0 0 flipped_loc]; %#ok<AGROW>
    end
end
out=flipud(out);
end
function state=DetermineSparseState(line)
% Determine the state of the sparser for every character to be able to detect comments and to split
% chars/strings and code.
% The input is a char vector, the output is a same-sized vector with values described below.
%
%state==-1 - escaping character for literal ' or "
%state==0  - normal code
%state==1  - opening a char array
%state==2  - inside a char array
%state==3  - closing a char array
%state==4  - opening a string
%state==5  - inside a string
%state==6  - closing a string

line=[' ' line ' '];%pad with space
idx1=strfind(line,'''');
idx2=strfind(line,'"');
state=zeros(size(line));

for n=sort([idx1 idx2])
    if state(n)==-1,continue,end % Escaping character for literal ' or ".
    if strcmp(line(n),'''')
        if state(n)==5,continue,end % Inside a string.
        if state(n)==2 % This either closes the char or starts a literal.
            if strcmp(line(n+1),'''')
                state(n+1)=-1;% Mark next character as escaping character.
            else
                state(n)=3;
                state((n+1):end)=0;
            end
        else
            % An apostrophe can only be a transpose if the previous character doesn't have anything
            % to do with a char array. So it is allowed to be a closing brace, a dot, a \w
            % character (i.e. a variable or function), a double quote (closing a string), or a
            % previous transpose.
            if state(n-1)~=3 && ~isempty(regexp(line(n-1),'[\]\)}.\w''"]')) %#ok<RGXP1>
                % The apostrophe is a transpose, so leave marked as normal code
            else
                % This starts a char array.
                state(n)=1;
                state((n+1):end)=2;
            end
        end
    else
        if state(n)==2,continue,end % Inside char array.
        if state(n)==5 % Either the closing character or the start of a literal.
            if strcmp(line(n+1),'"')
                state(n+1)=-1;% Mark next character as escaping character.
            else
                state(n)=6;
                state((n+1):end)=0;
            end
        else % This must be the start of a string.
            state(n)=4;
            state((n+1):end)=5;
        end
    end
end
state([1 end])=[];
end
function [SingleCommand,idx]=ExtractCommandsFromLine(str)
% Parse every line to single commands (the parts delimited by a comma or semicolon).
% This function assumes that the chars and strings are replaced by placeholders, or don't contain
% any brackets or braces.
temp=[' ' str(1:(end-1))];% Shift by 1 to make cumsum match the last element as well.
squareb = double(str=='[') - double(temp==']');
parenth = double(str=='(') - double(temp==')');
curlybr = double(str=='{') - double(temp=='}');
outer = ( cumsum(squareb)+cumsum(parenth)+cumsum(curlybr) ) == 0;
idx=find( outer & ( str==',' | str==';' ) );
if ~isempty(idx)
    % Gobble trailing spaces, commas, and semicolons.
    if all(ismember(str(idx(end):end),',; '))
        idx(end)=[];
    end
end
if isempty(idx)
    SingleCommand={str};
    idx=0;
else
    SingleCommand=cell(1+numel(idx),1);
    idx=[0 idx numel(str)+1];
    for m=1:numel(SingleCommand)
        SingleCommand{m}=str( (idx(m)+1):(idx(m+1)-1) );
    end
end
end
function pref=GetMaxLineLengthPref(varargin)
% Retrieve the maximum line length from the editor preferences. This isn't stored this in a
% persistent, in order to allow the user to change this between function calls.
%
% As far as I can tell there is no way to retrieve this programmatically for Octave, so the number
% mentioned in the style guide is used ( https://wiki.octave.org/Octave_style_guide ).
%
% This only works if the user has _ever_ opened the preference page. In general not opening the
% preference page will either result in the returned value being 0, or an empty value. For that
% reason a list of default values is included. This list may not be correct or complete, but it
% should be close enough, especially if the user doesn't care enough to have ever opened the
% preferences page.
persistent legacy default
if isempty(legacy)
    legacy.isOctave=exist('OCTAVE_VERSION', 'builtin');
    if ~legacy.isOctave
        legacy.settings=ifversion('>=','R2018a');
        legacy.com.mathworks=ifversion('>=','R13');
        legacy.R13=ifversion('==','R13');
    end
    if legacy.isOctave
        % https://wiki.octave.org/Octave_style_guide
        default=80;
    elseif ifversion('<=','R2010b')
        % The default for R2010b, R2010a and R14 seems to be 75.
        default=75;
    else
        % As far as I can tell the default is 80 at least since R2011a.
        default=80;
    end
end
if legacy.isOctave
    % As far as I can tell there is no way to retrieve this by code. This can be edited if a way is
    % found or is created.
    pref=0;
elseif legacy.settings
    % This was introduced in R2018a and should keep working after com.mathworks is removed.
    s = settings;pref=s.matlab.editor.displaysettings.linelimit.LineColumn.ActiveValue;
elseif legacy.com.mathworks
    % This works at least from R14SP3.
    pref=com.mathworks.services.Prefs.getIntegerPref('EditorRightTextLineLimit');%#ok<JAPIMATHWORKS>
elseif legacy.R13
    pref=com.mathworks.services.Prefs.getIntegerPref('EditorMaxCommentWidth');%#ok<JAPIMATHWORKS>
else % If you are using R12 or before, you have a lot more issues than just this.
    pref=80;
end
pref=double(pref);

if isempty(pref) || pref==0
    pref=default;
end
end
function tf=ifversion(test,Rxxxxab,Oct_flag,Oct_test,Oct_ver)
%Determine if the current version satisfies a version restriction
%
% To keep the function fast, no input checking is done. This function returns a NaN if a release
% name is used that is not in the dictionary.
%
% Syntax:
% tf=ifversion(test,Rxxxxab)
% tf=ifversion(test,Rxxxxab,'Octave',test_for_Octave,v_Octave)
%
% Output:
% tf       - If the current version satisfies the test this returns true.
%            This works similar to verLessThan.
%
% Inputs:
% Rxxxxab - Char array containing a release description (e.g. 'R13', 'R14SP2' or 'R2019a') or the
%           numeric version.
% test    - Char array containing a logical test. The interpretation of this is equivalent to
%           eval([current test Rxxxxab]). For examples, see below.
%
% Examples:
% ifversion('>=','R2009a') returns true when run on R2009a or later
% ifversion('<','R2016a') returns true when run on R2015b or older
% ifversion('==','R2018a') returns true only when run on R2018a
% ifversion('==',9.9) returns true only when run on R2020b
% ifversion('<',0,'Octave','>',0) returns true only on Octave
% ifversion('<',0,'Octave','>=',6) returns true only on Octave 6 and higher
%
% The conversion is based on a manual list and therefore needs to be updated manually, so it might
% not be complete. Although it should be possible to load the list from Wikipedia, this is not
% implemented.
%
%  _____________________________________________________________________________
% | Compatibility   | Windows XP/7/10 | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |-----------------|-----------------|------------------|----------------------|
% | ML R2020b       | W10: works      |  not tested      |  not tested          |
% | ML R2018a       | W10: works      |  works           |  not tested          |
% | ML R2015a       | W10: works      |  works           |  not tested          |
% | ML R2011a       | W10: works      |  works           |  not tested          |
% | ML R2010b       | not tested      |  works           |  not tested          |
% | ML R2010a       | W7:  works      |  not tested      |  not tested          |
% | ML 7.1 (R14SP3) | XP:  works      |  not tested      |  not tested          |
% | ML 6.5 (R13)    | W10: works      |  not tested      |  not tested          |
% | Octave 6.1.0    | W10: works      |  not tested      |  not tested          |
% | Octave 5.2.0    | W10: works      |  works           |  not tested          |
% | Octave 4.4.1    | W10: works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.0.5
% Date:    2020-12-08
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email = 'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

%The decimal of the version numbers are padded with a 0 to make sure v7.10 is larger than v7.9.
%This does mean that any numeric version input needs to be adapted. multiply by 100 and round to
%remove the potential for float rounding errors.
%Store in persistent for fast recall (don't use getpref, as that is slower than generating the
%variables and makes updating this function harder).
persistent  v_num v_dict octave
if isempty(v_num)
    %test if Octave is used instead of Matlab
    octave=exist('OCTAVE_VERSION', 'builtin');
    
    %get current version number
    v_num=version;
    ii=strfind(v_num,'.');if numel(ii)~=1,v_num(ii(2):end)='';ii=ii(1);end
    v_num=[str2double(v_num(1:(ii-1))) str2double(v_num((ii+1):end))];
    v_num=v_num(1)+v_num(2)/100;v_num=round(100*v_num);
    
    %get dictionary to use for ismember
    v_dict={...
        'R13' 605;'R13SP1' 605;'R13SP2' 605;'R14' 700;'R14SP1' 700;'R14SP2' 700;
        'R14SP3' 701;'R2006a' 702;'R2006b' 703;'R2007a' 704;'R2007b' 705;
        'R2008a' 706;'R2008b' 707;'R2009a' 708;'R2009b' 709;'R2010a' 710;
        'R2010b' 711;'R2011a' 712;'R2011b' 713;'R2012a' 714;'R2012b' 800;
        'R2013a' 801;'R2013b' 802;'R2014a' 803;'R2014b' 804;'R2015a' 805;
        'R2015b' 806;'R2016a' 900;'R2016b' 901;'R2017a' 902;'R2017b' 903;
        'R2018a' 904;'R2018b' 905;'R2019a' 906;'R2019b' 907;'R2020a' 908;
        'R2020b',909};
end

if octave
    if nargin==2
        warning('HJW:ifversion:NoOctaveTest',...
            ['No version test for Octave was provided.',char(10),...
            'This function might return an unexpected outcome.']) %#ok<CHARTEN>
        if isnumeric(Rxxxxab)
            v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
        else
            L=ismember(v_dict(:,1),Rxxxxab);
            if sum(L)~=1
                warning('HJW:ifversion:NotInDict',...
                    'The requested version is not in the hard-coded list.')
                tf=NaN;return
            else
                v=v_dict{L,2};
            end
        end
    elseif nargin==4
        % Undocumented shorthand syntax: skip the 'Octave' argument.
        [test,v]=deal(Oct_flag,Oct_test);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    else
        [test,v]=deal(Oct_test,Oct_ver);
        % Convert 4.1 to 401.
        v=0.1*v+0.9*fix(v);v=round(100*v);
    end
else
    % Convert R notation to numeric and convert 9.1 to 901.
    if isnumeric(Rxxxxab)
        v=0.1*Rxxxxab+0.9*fix(Rxxxxab);v=round(100*v);
    else
        L=ismember(v_dict(:,1),Rxxxxab);
        if sum(L)~=1
            warning('HJW:ifversion:NotInDict',...
                'The requested version is not in the hard-coded list.')
            tf=NaN;return
        else
            v=v_dict{L,2};
        end
    end
end
switch test
    case '==', tf= v_num == v;
    case '<' , tf= v_num <  v;
    case '<=', tf= v_num <= v;
    case '>' , tf= v_num >  v;
    case '>=', tf= v_num >= v;
end
end
function [names_in__fun,lines,position_in_line]=ListLocalFunctions(str)
%Extract the function names of all local functions and the line number they occur in.

isChar=isa(str,'char');if isChar,str={str};end

% Replace all string/char content with underscores to prevent incorrect matches.
[ignore,str]=SplitCodeAndChar(str); %#ok<ASGLU>

% Turn off the warning about the tokenize option being the default now.
w=warning('off','MATLAB:REGEXP:deprecated');

names_in__fun=cell(0);lines=zeros(0);position_in_line=zeros(0);n=0;
for line=1:numel(str)
    [SingleCommand,idx]=ExtractCommandsFromLine(str{line});
    for m=1:numel(SingleCommand)
        fun=ListLocalFunctions_helper(SingleCommand{m});
        if ~isempty(fun)
            n=n+1;
            names_in__fun{n,1}=fun;
            lines(n,1)=line;
            position_in_line(n,1)=idx(m);
        end
    end
end

%In a script this would fail, so assign a function name that is guaranteed not to cause problems
%later on: function. It can't be used as function name or variable name and will therefore not
%cause issues.
if numel(names_in__fun)==0,names_in__fun={'function'};end
if isChar,names_in__fun=names_in__fun{1};end
warning(w);% Reset warnings.
end
function fun=ListLocalFunctions_helper(str)
fun='';
RE='\s*function[\s\[](.*)';
ind=regexp(str,RE,'once');
if ~isempty(ind)
    temp=regexprep(str,RE,'$1', 'tokenize');
    % Since this is running on only a single command, we can take some shortcuts.    
    % Select everything past the equal sign, or from the start if there isn't any. Then we can crop
    % everything after the first parenthesis.
    temp=temp( max([1 1+strfind(temp,'=')]) : min([strfind(temp,'(')-1 numel(temp)]) );
    % Remove any stray spaces.
    fun=strrep(temp,' ','');
end
end
function var=ListVariables(str)
%List all variables that are defined in a cellstr.

% Turn off the warning about the tokenize option being the default now.
w=warning('off','MATLAB:REGEXP:deprecated');
persistent legacy
if isempty(legacy),try unique([1 1 2],'stable');legacy=false;catch,legacy=true;end,end

% Treat varargin and varargout as globals (i.e. do not return them as variable names).
var=cell(0);n=0;glo={'varargin';'varargout'};N=2;
for line=1:numel(str)
    SingleCommand=ExtractCommandsFromLine(str{line});
    
    for m=1:numel(SingleCommand)
        [var,n,glo,N]=ListVariables_ParseLine(SingleCommand{m},var,n,glo,N);
    end
end
if ~legacy
    var=unique(var,'stable');
else
    [ii,b]=unique(var); %#ok<ASGLU>
    var=var(sort(b));
end
var(ismember(var,glo))=[];% Remove globals from the list of variables.
warning(w)% Reset warning.
end
function [var,n,glo,N]=ListVariables_ParseLine(str,var,n,glo,N)
persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
temp=[' ' str];% Pad with space.

% Deal with function definition lines.
RE='\s*function[\s\[](.*)';
ind=regexp(str,RE,'once');
if ~isempty(ind)
    temp=regexprep(str,RE,'$1', 'tokenize');
    
    re='(^)|(\]?\s*=)\s*(([a-zA-Z]+[_a-zA-Z0-9]*))[^_a-zA-Z]';
    [a,b]=regexp(temp,re);
    try
        fun=temp(a(1):b(1));
        fun=regexprep(fun,'[^_a-zA-Z0-9]*([_a-zA-Z0-9]*)[^_a-zA-Z0-9]*','$1','tokenize');
    catch
        fun=regexprep(temp,re,'$1','tokenize');
    end
    temp=regexprep(temp,['\]?\s*=\s*' fun '\s*\((.*)\)'],',$1]', 'tokenize');
    
    temp=strrep(temp,fun,'');% This is needed for functions without inputs.
    if isempty(temp),return,end% No input and no output.
    [var,n]=ListVariables_ParseOutputVariables(temp,var,n);
    return
end

% Strip all bracketed characters, as they can't contain new variables.
L=zeros(size(temp));
for open ='{(',L(       temp==open  )= 1;end
for close='})',L([false temp==close])=-1;end
L=cumsum(L);temp(L(1:(end-1))>0)='';

% Check for persistent and global (but don't change the globals to avoid side-effects).
expr_={'[;,\s]' , '\s+(([a-zA-Z][a-zA-Z0-9_]*\s*)*)'};
expr= [expr_{1} 'persistent' expr_{2}];
if isOctave
    if isempty(regexp(temp,expr, 'once'))
        %No match, skip the parsing.
    else
        tokens=regexprep(temp,expr,'$1');
        [var,n]=ListVariables_ParseOutputVariables(tokens,var,n);
        temp=regexprep(temp,expr,'');
    end
else
    [s1,s2,tokens]=regexp(temp,expr);
    for k=numel(tokens):-1:1
        [var,n]=ListVariables_ParseOutputVariables(temp(tokens{k}(1):tokens{k}(2)),var,n);
        temp(s1(k):s2(k))='';
    end
end
expr= [expr_{1} 'global' expr_{2}];
if isOctave
    if isempty(regexp(temp,expr, 'once'))
        %No match, skip the parsing.
    else
        tokens=regexprep(temp,expr,'$1');
        [glo,N]=ListVariables_ParseOutputVariables(tokens,glo,N);
        temp=regexprep(temp,expr,'');
    end
else
    [s1,s2,tokens]=regexp(temp,expr);
    for k=numel(tokens):-1:1
        [glo,N]=ListVariables_ParseOutputVariables(temp(tokens{k}(1):tokens{k}(2)),glo,N);
        temp(s1(k):s2(k))='';
    end
end

% Replace comparators to make the rest of the parsing easier.
temp=strrep(temp,'==','<');
temp=strrep(temp,'~=','<');
temp=strrep(temp,'>=','<');
temp=strrep(temp,'<=','<');

for m1=numel(temp):-1:1
    % Find the last = and remove everything up to that point.
    ind=find(temp=='=');
    if isempty(ind)
        break
    else
        temp(ind(end):end)='';
        %hierarchy: [ , ; 1
        ind={max(find(temp=='[')), max([find(temp==',') find(temp==';')]), 1}; %#ok<MXFND>
        ind=[ind{:}];ind=ind(1);
        temp2=temp(ind(end):end);
        temp(ind(end):end)='';
        %temp2=regexprep(temp2,'\s+if\s+',''); % An if statement never leads to a new variable.
        temp2=regexprep(temp2,'\s+for\s+','');
        temp2=regexprep(temp2,'\s+try\s+','');
        [var,n]=ListVariables_ParseOutputVariables(temp2,var,n);
    end
end
end
function [var,n]=ListVariables_ParseOutputVariables(str,var,n)
% Strip outer square brackets.
[s1,s2,tokens]=regexp(str,'\[(.*)\]|((.*)?)'); %#ok<ASGLU>
str=str(tokens{1}(1):tokens{1}(2));

% Split variables.
expr = '[^a-zA-Z_\.]([a-zA-Z][a-zA-Z0-9_]*)';
% Pad with a space to allow filtering out of field names.
str=[' ' str];
[s1,s2,tokens]=regexp(str,expr); %#ok<ASGLU>
for k=1:numel(tokens)
    n=n+1;
    var{n,1}=str(tokens{k}(1):tokens{k}(2));
end
end
function str=MergeShortLines(str,max_length_after_merge,dont_merge_function_lines)
%Merge lines if they are short enough, adding padding spaces, commas or semicolons when needed.
if nargin<3,dont_merge_function_lines=true;end
pad=comma_or_semicolon(str);
for n=(numel(str)-1):-1:1
    line1=str{n};line2=str{n+1};
    % Strip leading whitespace.
    line2=regexprep(line2,'^\s*','');
    if dont_merge_function_lines
        % Function lines should not be merged.
        if strcmp(line2(1:(min(end,8))),'function'),continue,end
    end
    % Add a semicolon or comma to line1 to allow merging, but only if needed.
    if ~strcmp(line1(end),';') && ~strcmp(line1(end),',') 
        line1=[line1 pad(n)]; %#ok<AGROW>
    end
    % A try statement should be followed by a space or newline.
    if strcmp(line1((max(1,end-3)):max(1,end-1)),'try')
        line1(end)=' ';
    end
    % If the result is smaller than the maximum, merge the two lines.
    if numel(line1)+numel(line2) <= max_length_after_merge
        str{n}=[line1 line2];
        str(n+1)=[];
    end
end
end
function c=comma_or_semicolon(str)
% Return an array of spaces, commas and semicolons. A row contains a semicolon if the corresponding
% line ends inside a bracketed statement (either square or curly braces) and a comma otherwise. For
% the special case of an 'end' followed by a 'function', a space is returned.

% Use a comma by default.
c=repmat(',',numel(str),1);

% Find the position of the line ends in the merged string.
lengths=cumsum(cellfun('prodofsize',str));
str=horzcat(str{:});

% Use a semicolon if a line end occurs inside a bracketed statement.
squareb = double(str=='[') - double(str==']');
curlybr = double(str=='{') - double(str=='}');
outer = ( cumsum(squareb)+cumsum(curlybr) );
L=outer(lengths)~=0;
c(L)=';';

% Use a space (or a comma) to separate 'end' and 'function'. Just in case someone is using the
% Octave standard of 'endfunction', check if the location is a line end.
L=ismember(lengths,strfind(str,'endfunction')+2);
c(L)=SetDefault_end_function_char;
end
function varargout=regexp_outkeys(str,expression,varargin)
%Support the newer features of regexp in old releases (outkeys)
%
%Even if this extends functionality of regexp, cell array input is not supported in this function.
%
%Syntax:
% match=regexp_outkeys(str,expression,'match');
% split=regexp_outkeys(str,expression,'split');
% [match,split]=regexp_outkeys(str,expression,'match','split');
% [split,match]=regexp_outkeys(str,expression,'split','match');
%
%Example:
% str='lorem1 ipsum1.2 dolor3 sit amet 99 ';
% words=regexp_outkeys(str,' ','split')
% numbers=regexp_outkeys(str,'[0-9.]*','match')
%
%  _______________________________________________________________________
% | Compatibility | Windows 10  | Ubuntu 20.04 LTS | MacOS 10.15 Catalina |
% |---------------|-------------|------------------|----------------------|
% | ML R2020b     |  works      |  not tested      |  not tested          |
% | ML R2018a     |  works      |  works           |  not tested          |
% | ML R2015a     |  works      |  works           |  not tested          |
% | ML R2011a     |  works      |  works           |  not tested          |
% | ML 6.5 (R13)  |  works      |  not tested      |  not tested          |
% | Octave 5.2.0  |  works      |  works           |  not tested          |
% | Octave 4.4.1  |  works      |  not tested      |  works               |
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
%
% Version: 1.0.1 [not incremented, ver changed in table]
% Date:    2020-07-06
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email = 'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

if nargin<3
    error('HJW:regexp_outkeys:SyntaxError',...
        'No supported syntax used: at least 3 inputs expected.')
end
if ~(ischar(str) && ischar(expression))
    %extra params in varargin are checked inside the for-loop
    error('HJW:regexp_outkeys:InputError',...
        'All inputs must be char vectors.')
end

persistent legacy
if isempty(legacy)
    %The legacy struct contains the implemented options as field names. It is used in the error
    %message.
    %It is assumed that all Octave versions later than 4.0 support the expanded output, and all
    %earlier versions do not, even if it is very likely most versions will support it.
    
    %The switch to find matches was introduced in R14 (along with the 'tokenExtents', 'tokens' and
    %'names' output switches).
    legacy.match = ifversion('<','R14','Octave','<',4);
    
    %The split option was introduced in R2007b.
    legacy.split = ifversion('<','R2007b','Octave','<',4);
end
varargout=cell(size(varargin));%Pre-allocate output.
for param=1:(nargin-2)
    if ~ischar(varargin{param})
        error('HJW:regexp_outkeys:InputError',...
            'All inputs must be char vectors.')
    end
    switch lower(varargin{param})
        case 'match'
            if legacy.match
                %Legacy implementation.
                [s1,s2]=regexp(str,expression);
                match=cell(1,numel(s1));
                for n=1:numel(s1)
                    match{n}=str(s1(n):s2(n));
                end
            else
                match=regexp(str,expression,'match');
            end
            varargout{param}=match;
        case 'split'
            if legacy.split
                %Legacy implementation.
                [s1,s2]=regexp(str,expression);
                split=cell(1,numel(s1)+1);
                start_index=[s1 numel(str)+1];
                stop_index=[0 s2];
                for n=1:numel(start_index)
                    split{n}=str((stop_index(n)+1):(start_index(n)-1));
                end
            else
                split= regexp(str,expression,'split');
            end
            varargout{param}=split;
        otherwise
            fn=fieldnames(legacy);
            errorstr=['Extra regexp output type not implemented, only the following types are ',...
                'implemented:',char(10),sprintf('%s, ',fn{:})]; %#ok<CHARTEN>
            errorstr((end-1):end)='';%Remove trailing ', '
            error('HJW:regexp_outkeys:NotImplemented',errorstr)
    end
end
end
function end_function_char=SetDefault_end_function_char(end_function_char)
% ML6.5 requires a semicolon or comma between end and function, while newer release either require
% or allow a space. This function allows you to control which character is inserted.
% When called with an input, it will persistently store that input for future calls.
%
% Call with an empty input to set the default: a comma on ML6.5 and a space for the rest.
persistent ch
if nargin>0, ch=end_function_char;end
if isempty(ch)
    if ifversion('<',7,'Octave','>',0),ch=',';
    else,                              ch=' ';end
end
if nargout>0,end_function_char=ch;end
end
function functions=SortFunctions(str_,EntryFunctionNames)
%Sort the functions in the cellstr input alphabetically.
%The EntryFunctionNames will be put at the top and will not be sorted.
%If any code precedes a function definition on the same line, this will be split first.

[names__fun,lines,position_in_line]=ListLocalFunctions(str_);
if numel(lines)<=1 % Only one function (or none).
    functions={str_};return
end

% The position_in_line variable contains the starting indices of the command that contains the
% function definition. This index will be 1 even if there is whitespace before the function
% keyword, so any value above 1 means there is a command that should be moved to its own line.
if any(position_in_line>1)
    % Process all lines that need to be split
    inverted_list_of_lines=unique(lines(position_in_line>1)).';
    inverted_list_of_lines=inverted_list_of_lines(end:-1:1);
    for line=inverted_list_of_lines
        split_line=str_{line};
        pos=position_in_line(lines==line).';pos(pos==0)=[];
        lengths=diff([0 pos numel(split_line)]);
        split_line=mat2cell(split_line,1,lengths).';
        str_=[str_(1:(line-1));split_line;str_((line+1):end)];
    end
    
    % Re-run on the resulting cellstr. This is easier and less error-prone that correcting the
    % lines inside the previous loop.
    [names__fun,lines]=ListLocalFunctions(str_);
end    

% Find the order the functions should be sorted.
[entry_function_pos,order]=ismember(names__fun,EntryFunctionNames);
[ignore,order_of_non_entry_functions]=sort(names__fun(~entry_function_pos)); %#ok<ASGLU>
order(~entry_function_pos)=sum(entry_function_pos)+order_of_non_entry_functions;
[ignore,order]=sort(order); %#ok<ASGLU>

% Split functions to a cell array.
functions=cell(numel(lines),1);
lines=[lines;numel(str_)+1];
for n=1:numel(functions)
    functions{n}=str_(lines(n):(lines(n+1)-1));
end

% Apply sort.
functions=functions(order);
end
function [full_file,only_cmd]=SplitCodeAndChar(str)
% Split a cellstr into parts with code and parts with char/string definitions.
% The first output is an array with the even columns containing the code that creates chars/strings
% and the odd columns the rest of the code.
% The second output is a column vector cell where all char/string content is replaced with
% underscores, so the syntax is easier to parse.
persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
full_file=cell(numel(str),1);
only_cmd=cell(numel(str),1);
for n=1:numel(str)
    % Split each line into code and char/string.
    temp=SplitLineToCodeAndChar(str{n});
    % Store the code in the odd columns and the chars/strings in the even columns.
    full_file(n,1:2:(2*size(temp,2)))=temp(1,:);
    full_file(n,2:2:(2*size(temp,2)-1))=temp(2,1:(end-1));
    for m=1:(size(temp,2)-1)
        % Replace with dummy arrays to simplify syntax handling.
        quot=temp{2,m}(1);%This is either a single or a double quote.
        temp{2,m}=[quot repmat('_',1,numel(temp{2,m})-2) quot];
    end
    if isOctave,temp(cellfun('isempty',temp))={''};end
    temp=[temp{:}];
    only_cmd{n}=temp;
end
end
function [ostr,CRLF,hasTrailingLF]=SplitLines(istr)
%Split a char array or string scalar to a cellstr.
istr=char(istr);
hasTrailingLF=false;
if ~( any(istr==10) || any(istr==13) )
    if ispc,CRLF=char([13 10]);else,CRLF=char(10);end %#ok<CHARTEN>
    ostr={istr};
else
    % Find the newline style: \n, \r\n, \n\r, or \r.
    CRLF_ind=find(istr==10 | istr==13);CRLF_ind=CRLF_ind(1:min(2,end));
    if numel(CRLF_ind)==2 && ...
            ( diff(CRLF_ind)~=1 || istr(CRLF_ind(1))==istr(CRLF_ind(2)) )
        % Revert to \n or \r.
        CRLF_ind(2)=[];
    end
    CRLF=istr(CRLF_ind);
    % Determine if the char contains a trailing newline.
    if numel(istr)>1
        hasTrailingLF=strcmp(CRLF,istr((end-numel(CRLF)+1):end));
    end
    
    % Split over the line ends.
    ostr=regexp_outkeys(istr,CRLF,'split');
    ostr=ostr(:);
end
end
function split_line=SplitLineToCodeAndChar(str)
%Split a char array into code and char components.
%Syntax:
%  split_line=SplitLineToCodeAndChar(str)
%
% str:        1xn char array with a single line of Matlab code
% split_line: 2xm cell array with all chars on the second row (such that str==horzcat(full_line{:})

str=[' ' str ' '];% Pad with spaces to standardize state detection output shape.

% Determine which characters are inside chars/strings.
state=DetermineSparseState(str);

y=diff(state==0);
CodeStarts=[1 find(y==1)+1             ];
CodeEndsAt=[  find(y==-1)  numel(state)];
split_line=cell(2,numel(CodeStarts));
for n=1:numel(CodeStarts)
    % Fill the code part.
    a=CodeStarts(n);b=CodeEndsAt(n);
    split_line{1,n}=str(a:b);
    % Fill the char/string part.
    a=CodeEndsAt(n)+1;b=CodeStarts(min(n+1,end))-1;
    split_line{2,n}=str(a:b);% This will be empty if a>b.
end

% Remove padding.
split_line{1, 1 }( 1 )='';
split_line{1,end}(end)='';
end
function ostr = StripComments(istr)
%Strip comments from a char array or string with MATLAB code.
%
%Syntax:
% ostr = StripComments(istr)
%
% istr can be either a char array or a cellstr.
% ostr will have the same class as istr.
%
% This function will strip all comments:
%  - all lines from %{ up to and including %}
%  - all text after ... (a space will be added only if needed for a valid syntax)
%  - all text after % (including % itself)
% All trailing whitespace will also be removed.

OutputString=isa(istr,'string');
if OutputString,istr=cellstr(istr);end
if isa(istr,'cell')
    for n=1:numel(istr)
        % All character arrays must be row vectors.
        if (ndims(istr{n}) > 2) || (any(size(istr{n}) > 0) && (size(istr{n}, 1) ~= 1)) %#ok<ISMAT>
            error('HJW:StripComments:InvalidInput','All character arrays must be row vectors.');
        end
    end
end

% If input is a character array, put it into a cell array.  We'll later make
% sure output is a character array if input is a character array.
OutputChar=ischar(istr);

% Make sure istr is a proper cellstr. SplitLines will split a char into a cellstr.
if OutputChar
    [istr,CRLF,hasTrailingLF]=SplitLines(istr);
else
    CRLF='';hasTrailingLF=false;% Set default.
    for n=1:numel(istr)
        [istr{n},CRLF_,hasTrailingLF_]=SplitLines(istr{n});
        hasTrailingLF=hasTrailingLF || hasTrailingLF_;
        if numel(istr{n})>1,CRLF=CRLF_;end
    end
    if isempty(CRLF),CRLF=CRLF_;end % Use default if never set.
    % Re-concatenate the output.
    istr=vertcat(istr{:});
end

% Remove block comments.
block_comment=cell(1,2);
for n=1:2
    if n==1,brace='\{';else,brace='\}';end
    % Block comments only allow whitespace on the lines with %{ and %}.
    block_comment{n}=regexp(istr,['^[^\S]*%' brace '[^\S]*$']);
    if isa(block_comment{n},'double'),block_comment{n}=block_comment(n);end % Only needed pre-R14.
    block_comment{n}=~cellfun('isempty',block_comment{n});
end
% Block comments can be nested. The cumsum skips the last closing %}.
block_comment=cumsum(block_comment{1}-block_comment{2})>0 | block_comment{2};
% Remove block comments from the cell array.
istr(block_comment)=[];

% Initialize output.
ostr = istr;
% Iterate over each element in the cell array.
for n = 1:numel(ostr)
    % Get the ith input string.
    str = ostr{n};
    
    % Remove a comment on a line possibly containing MATLAB code.
    str = StripComments__line(str);
    
    % Remove trailing whitespace.
    IsWhitespace=isspace(str);if all(IsWhitespace),idx=0;else,idx=find(~IsWhitespace);end
    if ~isempty(idx)&&idx(end)<numel(str),str((idx(end)+1):end)=[];end
        
    % Store in the output cell array.
    ostr{n} = str;
end

% Loop through the line to remove line continuations (back to front so lines can be deleted).
for n = numel(ostr):-1:1
    removed=istr{n}((1+numel(ostr{n})):end);
    % Does the removed text start with ... (optionally padded with whitespace)?
    [ind_s,ind_e]=regexp(removed,'[\s]*[^\s\.%]?[\s]*\.{3}');
    if ~isempty(ind_e)&&ind_s(1)==1&&numel(istr)>n
        % Find the index to strip leading whitespace from the next line.
        idx=find(~isspace(ostr{n+1}));
        if ~isempty(idx)&&idx(1)>1,idx=idx(1);else,idx=1;end
        % Merge, adding a space (only if required) to avoid syntax errors.
        if numel(ostr{n})>=1 && any(ostr{n}(end)=='~,;+-*/^\@<>&|=')
            padding='';else,padding=' ';
        end
        if idx==1 % Avoid indexing if we don't need to.
            ostr{n}=[ostr{n} padding ostr{n+1}];ostr(n+1)=[];
        else
            ostr{n}=[ostr{n} padding ostr{n+1}(idx:end)];ostr(n+1)=[];
        end
    end
end
% Remove blank lines.
ostr(cellfun('isempty',ostr))=[];
% If the input was a character array, make sure output is so too.
if OutputChar
    ostr = sprintf(['%s' CRLF],ostr{:});
    if ~hasTrailingLF
        try ostr((end-numel(CRLF)+1):end) = '';catch,end
    end
elseif OutputString
    ostr=string(ostr);
end
end
function stripped=StripComments__line(line)
%Strip comments from a char array or cellstr with MATLAB code.

[str,CRLF]=SplitLines(line);
if numel(str)>1
    % Parse each line separately.
    for n=1:numel(str)
        str{n}=StripComments__line(str{n});
    end
    % Add the line end and make the last cell a char to prevent a data type error.
    str=str.';str(2,1:(end-1))={CRLF};str{end}='';
    % Merge back to a single char.
    stripped=horzcat(str{:});
    % Cast to string if the input was a string
    if isa(line,'string'),stripped=string(stripped);end
    return
else
    line=str{1};
end

% Determine which characters are inside chars/strings.
state=DetermineSparseState(line);

% Find the first percent symbol or ellipsis that has a state of 0 (i.e. code).
ind1=strfind(line,'%');ind2=strfind(line,'...');
ind=[ind1(state(ind1)==0) ind2(state(ind2)==0)];
if ~isempty(ind),line(min(ind):end)='';end
stripped=line;
end
function [isLogical,val]=test_if_scalar_logical(val)
%Test if the input is a scalar logical or convertible to it.
%(use the first output to trigger an input error, use the second as the parsed input)
%
% Allowed values:
%- true or false
%- 1 or 0
%- 'on' or 'off'
%- matlab.lang.OnOffSwitchState.on or matlab.lang.OnOffSwitchState.off
persistent states
if isempty(states)
    states={true,false;...
        1,0;...
        'on','off'};
    try
        states(end+1,:)=eval('{"on","off"}');
    catch
    end
end
isLogical=true;
try
    for n=1:size(states,1)
        for m=1:2
            if isequal(val,states{n,m})
                val=states{1,m};return
            end
        end
    end
    if isa(val,'matlab.lang.OnOffSwitchState')
        val=logical(val);return
    end
catch
end
isLogical=false;
end