function aaa___minify___test(varargin)
% Test the most important stages separately, and check if the output of the whole chain is stable.

global HJW___test_suit___debug_hook_data
HJW___test_suit___debug_hook_data=[];% Reset so it doesn't interfere with the first runs.

%Test the comment stripping.
stage='comment stripping';
[str,ground_truth]=comment_stripper_test_data;
try ME=[]; %#ok<NASGU>
    stripped = StripComments(str);
    failflag=false;
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    disp(ME.message),disp(ME.identifier)
    failflag=true;
end
if failflag || ( ~isequal(stripped,ground_truth)  )%&& ifversion('>',7.1,'Octave','>',0)
    error(['test did not match expected output (' stage ')'])
end

%Test the char/string detection.
stage='char/string detection';
[str,ground_truth]=char_and_code_split_test_data;
try ME=[]; %#ok<NASGU>
    full_file=SplitCodeAndChar(str);
    
    % A 0x0 double and a 1x0 char should not be treated differently, as the result is equivalent.
    ground_truth(cellfun('isempty',ground_truth))={''};
    full_file(   cellfun('isempty',full_file   ))={''};
    failflag=false;
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    disp(ME.message),disp(ME.identifier)
    failflag=true;
end
if failflag || ~isequal(full_file,ground_truth)
    error(['test did not match expected output (' stage ')'])
end

%Test the function and variable name detection.
stage='function/variable name detection';
[str,var_ground_truth,fun_ground_truth]=generate_fill_workspace;
try ME=[]; %#ok<NASGU>
    var=ListVariables(str);
    fun=ListLocalFunctions(str);
    failflag=false;
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    disp(ME.message),disp(ME.identifier)
    failflag=true;
end
if failflag || (~isequal(var,var_ground_truth) || ~isequal(fun,fun_ground_truth) )
    error(['test did not match expected output (' stage ')'])
end

%Test the whole chain. The minified code should be stable for a particular version of the code and
%the specific settings that are chosen. A newer version might return a different minified result,
%and therefore a different checksum.
str=[generate_fill_workspace;char_and_code_split_test_data;comment_stripper_test_data];
FormatSpec='omnibus test part %d did not match expected hash\n(expected %s, returned %s)';

part=1;ME=[];H='';H_=''; %#ok<NASGU>
try
    str_=minify(str,true,80,{'fill_workspace','fun_loc_1'},' ');
    H=ComputeNonCryptHash(str_,128,'-v2');
    H_='7E9122EF73B68A13258675C9C892F679';
    failflag=false;
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    disp(ME.message),disp(ME.identifier)
    failflag=true;
end
if failflag || ~strcmp(H,H_)
    error(FormatSpec,part,H_,H)
end

part=2;ME=[];H='';H_=''; %#ok<NASGU>
BeastMode=struct('compress_to_block',true);
try
    str_=minify(str,BeastMode,80,{'fill_workspace','fun_loc_1'},' ');
    H=ComputeNonCryptHash(str_,128,'-v2');
    H_='FC671B56444BC0D4764B3B77911DA697';
    failflag=false;
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    disp(ME.message),disp(ME.identifier)
    failflag=true;
end
if failflag || ~strcmp(H,H_)
    error(FormatSpec,part,H_,H)
end

part=3;ME=[];H='';H_=''; %#ok<NASGU>
BeastMode=true;
% Use the debug hook to set every local function as an entry function.
HJW___test_suit___debug_hook_data=struct('action','return','data',{{true}});
try
    str_=minify(str,BeastMode,80,[],' ');
    H=ComputeNonCryptHash(str_,128,'-v2');
    H_='53AA0EC79DFD324324C2B5A0D2676EE8';
    failflag=false;
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    disp(ME.message),disp(ME.identifier)
    failflag=true;
end
if failflag || ~strcmp(H,H_)
    error(FormatSpec,part,H_,H)
end

part=4;ME=[];H='';H_=''; %#ok<NASGU>
try
    str_=minify(str,true,80,{'fill_workspace'},' ');
    H=ComputeNonCryptHash(str_,128,'-v2');
    H_='E362A42F3F8EC44A5742C3CACDC3B629';
    failflag=false;
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    disp(ME.message),disp(ME.identifier)
    failflag=true;
end
if failflag || ~strcmp(H,H_)
    error(FormatSpec,part,H_,H)
end

part=5;ME=[];H='';H_=''; %#ok<NASGU>
try
    str_=minify(str,false,80,{'fill_workspace'},' ');
    H=ComputeNonCryptHash(str_,128,'-v2');
    H_='B985E20AFB7B5FB62630EE0561BF254F';
    failflag=false;
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    disp(ME.message),disp(ME.identifier)
    failflag=true;
end
if failflag || ~strcmp(H,H_)
    error(FormatSpec,part,H_,H)
end

part=6;ME=[];H='';H_=''; %#ok<NASGU>
BeastMode=struct('keep_original_function_names',true);
try
    str_=minify(str,BeastMode,80,{'fill_workspace'},' ');
    H=ComputeNonCryptHash(str_,128,'-v2');
    H_='53AA0EC79DFD324324C2B5A0D2676EE8';
    failflag=false;
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    disp(ME.message),disp(ME.identifier)
    failflag=true;
end
if failflag || ~strcmp(H,H_)
    error(FormatSpec,part,H_,H)
end

part=7;ME=[];H='';H_=''; %#ok<NASGU>
BeastMode=struct('contains_nested_functions',true);
try
    str_=minify(str,BeastMode,80,{'fill_workspace'},' ');
    H=ComputeNonCryptHash(str_,128,'-v2');
    H_='3788C14F8A53F72CD5C02BE5767F2F8C';
    failflag=false;
catch ME;if isempty(ME),ME=lasterror;end %#ok<LERR>
    disp(ME.message),disp(ME.identifier)
    failflag=true;
end
if failflag || ~strcmp(H,H_)
    error(FormatSpec,part,H_,H)
end


disp('test finished successfully')
end
function tf=transposing_scalar_string_is_fixed
% Transposing a scalar string (without a dot) did not work with the original version. This function
% controls the addition of a few test cases. Even though this code should not occur in normal code
% (as it has no effect), the comment stripper was rewritten from scratch, which also resolved this
% issue.
persistent tf_
if isempty(tf_)
    str='""''';
    tf_=isequal(StripComments(str),str);
    if ~tf_,error('This should have been fixed already.'),end
end
tf=tf_;
end
function [str,var_ground_truth,fun_ground_truth]=generate_fill_workspace
%Line continuations (i.e. ...) and comments should not be tested with this, as they should be have
%been removed by the time the function and variable name detection is performed.
str={...
    '%This file lists most possible syntaxes that result in workspace variables.';...
    '%All variables named var__ should be found, but none of the old_var__ variables.';...
    'function [var01,var02]=fill_workspace(var03,var04,var05)';...
    ' function varargout=fun_loc_1(varargin);';...
    '[var06,var07]=fun_ext_1(old_var1)';...
    'var08=fun_ext_2(old_var2);[var09,var10]=fun_ext_3(old_var3);';...
    'var11.field1=fun_ext_4;';...
    'var12.(old_var4)=fun_ext_5;';...
    'persistent var13;var14=fun_ext_6,persistent var15;';...
    ';persistent var16 var17';...
    'global var_glo_01';...
    'persistent var18;[old_var5];global var_glo_02;';...
    'old_var6==old_var7';...
    '[old_var8 old_var9]';...
    'if old_var10';...
    '    1+1;';...
    'end';...
    '  for var19=1:numel(old_var_11),end,if fun_ext_7(old_var11),var20=fun_ext_8(old_var12);end';...
    'if fun_ext_9(old_var13)';...
    'end';...
    'old_var14,function fun_loc_2';...
    'function[var21]=fun_loc_3';...
    'function var22=fun_loc_4';...
    'function [var23,var24]=fun_loc_5,fun_ext_10(old_var15)';...
    'persistent persistent_var_25';...
    'persistent var_26_persistent';...
    'persistent var_persistent_27';...
    'persistent_var28=old_var_16;';...
    'function [var29,var30]=fun_loc_6,[old_var17(old_var18)];function var31=fun_loc_7';...
    '@fun_loc_6'};
var_ground_truth={'var01';'var02';'var03';'var04';'var05';'var06';'var07';'var08';'var09';...
    'var10';'var11';'var12';'var13';'var14';'var15';'var16';'var17';'var18';'var19';'var20';...
    'var21';'var22';'var23';'var24';'persistent_var_25';'var_26_persistent';...
    'var_persistent_27';'persistent_var28';'var29';'var30';'var31'};
fun_ground_truth={'fill_workspace';'fun_loc_1';'fun_loc_2';'fun_loc_3';'fun_loc_4';'fun_loc_5';...
    'fun_loc_6';'fun_loc_7'};
end
function [str,ground_truth]=char_and_code_split_test_data
str={'data=[bit1'';bit2'';bit3'';bit4''];';
    'data=[bit1'';''foo%'';bit2'';bit3'';bit4''];';
    'data=[bit1'';"foo%";bit2'';bit3'';bit4''];';
    'data=[bit1'';''foo'';''bar'';bit2.'';bit3'';bit4''];';
    '''foo'';''bar'';bit2.'';bit3'';[bit4''];';
    ' ''foo'';''bar'';bit2.'';bit3'';[bit4''];';
    'data=[bit1'';"foo%".'';bit2'';bit3'';bit4''];';
    'data=[bit1'';"foo%"'';bit2'';bit3'';bit4''];'};
ground_truth={'data=[bit1'';bit2'';bit3'';bit4''];',[],[],[],[];...
    'data=[bit1'';','''foo%''',';bit2'';bit3'';bit4''];',[],[];...
    'data=[bit1'';','"foo%"',';bit2'';bit3'';bit4''];',[],[];...
    'data=[bit1'';','''foo''',';','''bar''',';bit2.'';bit3'';bit4''];';...
    [],'''foo''',';','''bar''',';bit2.'';bit3'';[bit4''];';...
    ' ','''foo''',';','''bar''',';bit2.'';bit3'';[bit4''];';
    'data=[bit1'';','"foo%"','.'';bit2'';bit3'';bit4''];',[],[];
    'data=[bit1'';','"foo%"',''';bit2'';bit3'';bit4''];',[],[]};
if ~transposing_scalar_string_is_fixed
    str(8)=[];
    ground_truth(8,:)=[];
end
end
function [str,ground_truth]=comment_stripper_test_data
str={'%just a few examples of hard to parse syntaxes';...
    '''%".''''...''%foo';...
    '%{';...
    'bar';...
    '%}';...
    '"foo"';...
    'bar.''''';...
    '[foo...';...
    '    bar]';...
    '[foo...foobar';...
    '    bar]';...
    '""""''%''...';...
    'foo';...
    '[""]''.''%foo';...
    '[1 2';...
    '    3 4]';...
    '%{foo}';...
    'bar'
    '   %foo'};
ground_truth={...
    '''%".''''...''';...
    '"foo"';...
    'bar.''''';...
    '[foo bar]';...
    '[foo bar]';...
    '""""''';...
    'foo';...
    '[""]''.''';...
    '[1 2';...
    '    3 4]';...
    'bar'};
if ~transposing_scalar_string_is_fixed
    str([12 13])=[];
    ground_truth(6)=[];
end
end
function out=bsxfun_plus(in1,in2)
%Implicit expansion for plus(), but without any input validation.
try
    out=in1+in2;
catch
    try
        out=bsxfun(@plus,in1,in2);
    catch
        sz1=size(in1);                    sz2=size(in2);
        in1=repmat(in1,max(1,sz2./sz1));  in2=repmat(in2,max(1,sz1./sz2));
        out=in1+in2;
    end
end
end
function data=cast_to_uint16_vector(data,re_encode_char,string_to_cellstr)
%Linearize the input data and convert it to a uint16 vector.
if isa(data,'uint16')
    %Append the array size and type to make it influence the hash.
    c='uint16';sz=size(data).';
    data=reshape(data,[],1);%linearize
    data=[data;uint16(c.');uint16(mod(sz,2^16))];
    return
end
data=cast_to_uint16_vector__cell({data},re_encode_char,string_to_cellstr);
data([end-1 end])=[];%Remove the [1 1] that is added because of the wrapping in a cell.
end
function data=cast_to_uint16_vector__cell(data,re_encode_char,string_to_cellstr)
sz=size(data).';data=data(:);
for n=1:numel(data)
    if numel(data{n})==0
        c=double(class(data{n})');
        data{n}=uint16([0;c;size(data{n})']);
        continue
    end
    switch class(data{n})
        case {'double','single'}
            data{n}=cast_to_uint16_vector__floats(data{n});
        case 'logical'
            data{n}=cast_to_uint16_vector__logical(data{n});
        case {'uint8','uint16','uint32','uint64','int8','int16','int32','int64'}
            data{n}=cast_to_uint16_vector__integer(data{n});
        case 'char'
            data{n}=cast_to_uint16_vector__char(data{n},re_encode_char);
        case 'string'
            data{n}=cast_to_uint16_vector__string(data{n},re_encode_char,string_to_cellstr);
        case 'cell'
            data{n}=cast_to_uint16_vector__cell(data{n},re_encode_char,string_to_cellstr);
        case 'struct'
            data{n}=cast_to_uint16_vector__struct(data{n},re_encode_char,string_to_cellstr);
        otherwise
            error('HJW:cast_to_uint16_vector:nosupport',...
                'Unsupported data type in nested variable')
    end
end
data=cell2mat(data);%Merge all cell contents.
data=[data;uint16(mod(sz,2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__floats(data)
sz=size(data).';c=class(data);%The rest of the function treats singles as double.

%Convert to a uint64, separate it into 4 words and everything merge into a vector.
[bit64,bit4]=typecast_double_uint64(double(data));
bit4_round=mod(bit64,2^16);bit64=bit64-bit4_round;bit64=bit64/2^16;bit4=bit4.';
bit3      =mod(bit64,2^16);bit64=bit64-bit3;      bit64=bit64/2^16;bit3=bit3.';
bit2      =mod(bit64,2^16);bit64=bit64-bit2;      bit64=bit64/2^16;bit2=bit2.';
bit1      =mod(bit64,2^16);                                        bit1=bit1.';
data=[bit1;bit2;bit3;bit4];
data=uint16(data(:));

%Append the array size and type to make it influence the hash.
data=[data;uint16(c.');uint16(mod(sz,2^16))];
end
function data=cast_to_uint16_vector__logical(data)
sz=size(data).';data=data(:);

if mod(numel(data),16) %Pad to 16 bits with 0.
    data(16*ceil(numel(data)/16))=0;
end
vector=uint16(2.^(15:-1:0))';
data=uint16(reshape(data,16,[]));
try
    data=data.*vector;
catch %No implicit expansion.
    %times() is not defined for uint16 on ML6.5
    data=double(data).*repmat(double(vector),[1 size(data,2)]);
    data=uint16(data);
end
data=uint16(sum(data,1)).';

data=[data;uint16(mod(sz,2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__integer(data)
%Large values (>2^52) will not have integer precision due to a conversion to double.
%This conversion is done, because of limited availability of operations on ML6.5.
sz=size(data).';data=data(:);

c=class(data);
if c(1)~='u'
    %Shift int* values to the uint* data range
    data=double(data)-double(eval([c '(-inf)']));
else
    data=double(data);
end
switch c(end)
    case '8'
        %Append a 0 for odd length and merge pairs.
        if mod(numel(data),2),data(end+1)=0;end
        data=reshape(data,[],2);
        data=data(:,1)*255+data(:,2);
        data=uint16(data);
    case '6'
        data=uint16(data);
    case '2'
        %Split to 2 words.
        bit1=floor(data/2^16);bit1=bit1.';
        bit2=mod(data,2^16);  bit2=bit2.';
        data=[bit1;bit2];
        data=uint16(data(:));
    case '4'
        %Split to 4 words. Note that bit4 contains a rounding error for data >2^52.
        bit64=data;
        bit4_round=mod(bit64,2^16);bit64=bit64-bit4_round;bit64=bit64/2^16;bit4_round=bit4_round.';
        bit3      =mod(bit64,2^16);bit64=bit64-bit3;      bit64=bit64/2^16;bit3=bit3.';
        bit2      =mod(bit64,2^16);bit64=bit64-bit2;      bit64=bit64/2^16;bit2=bit2.';
        bit1      =mod(bit64,2^16);                                        bit1=bit1.';
        
        data=[bit1;bit2;bit3;bit4_round];
        data=uint16(data(:));
end
%Append the array size and type to make it influence the hash.
data=[data;uint16(c.');uint16(mod(sz,2^16))];
end
function data=cast_to_uint16_vector__char(data,re_encode_char)
persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
if isOctave && re_encode_char
    isColVector = size(data,1)==numel(data);
    if isColVector,data=data.';end
    data=cellstr(data);
    % Decode from UTF-8 and encode with UTF-16 (the output will be uint16).
    for n=1:numel(data)
        data{n}=unicode_to_char(UTF8_to_unicode(data{n},true));
    end
    % Pad with spaces (but in uint16).
    cell_length=cellfun('length',data);longest_cell=max(cell_length);
    for n=find(cell_length<longest_cell)
        data{n}( (numel(data{n})+1) : longest_cell)=uint16(' ');
    end
    data=cell2mat(data);
    if isColVector,data=data.';end
end
sz=size(data).';data=data(:);
data=uint16(data);%Chars are 16 bit in Matlab, as they are encoded with UTF-16.
data=[data;uint16(mod(sz,2^16))];%Append the array size to make it influence the hash.
end
function data=cast_to_uint16_vector__string(data,re_encode_char,string_to_cellstr)
if string_to_cellstr
    data=cellstr(data);
    data=cast_to_uint16_vector__cell(data,re_encode_char,string_to_cellstr);
else
    data=char(data);%Cast to char instead of a cell array of chars.
    data=cast_to_uint16_vector__char(data,re_encode_char);
end
end
function data=cast_to_uint16_vector__struct(data,re_encode_char,string_to_cellstr)
sz=size(data).';data=data(:);
fn=fieldnames(data);
output=cell(2,numel(fn));
for n=1:numel(fn)
    output{1,n}=fn{n};
    output{2,n}={data.(fn{n})};
end
data=cast_to_uint16_vector__cell(output,re_encode_char,string_to_cellstr);
data=[data;uint16(mod(sz,2^16))];%Append the array size to make it influence the hash.
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
function hash=ComputeNonCryptHash(data,varargin)
%Compute a non-cryptographic hash
%
% This function is intended to be fast, but without requiring a Java or mex implementation to do
% the actual hashing. It was *not* checked for any security flaws and is therefore probably
% vulnerable to most attacks.
% Non-cryptographic hashes should only be used as a checksum. Don't use this to do things like
% storing passwords.
%
%syntax:
%  hash=ComputeNonCryptHash(data)
%  hash=ComputeNonCryptHash(data,HashLength)
%  hash=ComputeNonCryptHash(___,VersionFlag)
%
%data        The data to be hashed. Most common data types are allowed: uint*, int*, char, cell,
%            struct, double, or single (string is cast to a cell array of chars). The contents of
%            the nested data types (i.e. cell and struct) must also be one of the mentioned data
%            types.
%HashLength  The length of the hash (the number of bits). This value must be a multiple of 16. The
%            default is 256 bits. Depending on your input 64 bits might have some collisions, but
%            64 bits and higher should be safe.
%VersionFlag Either '-v1', '-v2'. This is provided for backwards compatibility. Version 1 of this
%            function has many hash collisions for scalar doubles and attempts to cast strings to
%            chars, instead of casting to a cell array of chars. Version 2 also decodes the UTF-8
%            chars from Octave and re-encodes them with UTF-16. That way the output is stable for
%            the Unicode code points.
%
%hash        The hash in an upper case hexadecimal char vector of size 1x(HashLength/4).
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
% Version: 2.0.1
% Date:    2020-12-09
% Author:  H.J. Wisselink
% Licence: CC by-nc-sa 4.0 ( https://creativecommons.org/licenses/by-nc-sa/4.0 )
% Email = 'h_j_wisselink*alumnus_utwente_nl';
% Real_email = regexprep(Email,{'*','_'},{'@','.'})

if nargin<1
    error('HJW:ComputeNonCryptHash:InputIncorrect','At least 1 input required.')
end
if nargin<2
    HashLength=256;
    UseVersion=2;
elseif nargin<3
    HashLength=varargin{1};
    UseVersion=2;
else
    HashLength=varargin{1};
    UseVersion=varargin{2};
    try
        if isa(UseVersion,'string'),UseVersion=char(UseVersion);end
        if ~isa(UseVersion,'char'), error('trigger'); end
        UseVersion=str2double(UseVersion(3:end));
        if isnan(UseVersion) || round(UseVersion)~=UseVersion || UseVersion>2
            error('trigger');
        end
    catch
            error('HJW:ComputeNonCryptHash:InputIncorrect',...
                'Version input incorrect. Must be ''-v1'', ''-v2''.')
    end
end
if numel(HashLength)~=1 || ~isnumeric(HashLength) || mod(HashLength,16)~=0 || HashLength<16
    error('HJW:ComputeNonCryptHash:InputIncorrect',...
        'Second input (hash length) must be a multiple of 16.')
end
    
try
    %Convert the input to an uint16 array (Nx1).
    re_encode_char_on_Octave=UseVersion>=2;%Chars will be decoded and encoded with UTF-16.
    string_to_cellstr=UseVersion>=2;%Cast strings to a cell of chars, instead of a char.
    data=cast_to_uint16_vector(data,re_encode_char_on_Octave,string_to_cellstr);
catch
    ME=lasterror; %#ok<LERR>
    if strcmp(ME.identifier,'MATLAB:nomem')
        %rethrow memory error
        rethrow(ME)
    else
        error('HJW:ComputeNonCryptHash:UnwindFailed',...
            'The nested input contains an unsupported data type.')
    end
end

%Extend to a multiple of HashLength bits. Padding with zeros is generally not advised, and the
%performance penalty for this extension (compared to padding with zeros) should be negligible.
if mod(numel(data),HashLength/16)
    extra=uint16(1:HashLength/16).'; extra(1:mod(numel(data),HashLength/16))=[];
    data=[data;extra];
end

%Add perturbation to the data and convert to 16xN logical. Then further perturb the intermediate
%result by circshifting every column (by col_index-1 positions) and doing an XOR in blocks with
%itself (by reshaping and transposing).
if UseVersion==1
    data=ComputeNonCryptHash_shuffle_uint16(data);
    data=ComputeNonCryptHash_uint16_to_logical(data);
    data=xor(data,reshape(data,[],16).');
else
    data=ComputeNonCryptHash_shuffle_uint16(data);
    data=ComputeNonCryptHash_uint16_to_logical(data);
    data=circshift_by_col(data);
    %data=circshift_by_col(data.').';
    %data=xor(data,not(reshape(data,[],16).'));
    %data=checker_shift(data);
end

%Reshape to HashLength cols and collapse the key size down to the hash length by counting the
%number of true bits (even=1, odd=0).
data=mod(sum(reshape(data,HashLength,[]),2),2);
data=ComputeNonCryptHash_logical_to_uint16(data);

if nargin>3
    hash=data;%Return uint16 for the salting.
    return
end

%Perturb the hash, analogous to salting. This function computes the hash of the hash and applies a
%few operations to the data to increase the randomness of the end result.
data=ComputeNonCryptHash_add_salt(data,UseVersion);

%Convert the (HashLength/16)x1 uint16 to a hash string by encoding it as hexadecimal.
hash=ComputeNonCryptHash_dec2hex(data);hash=reshape(hash.',1,[]);
end
function data=circshift_by_col(data)
%Circshift every column by col_index-1 positions.
persistent LUT
sz=size(data);
if isempty(LUT) || any(size(LUT)<sz) || isempty(LUT{sz(1),sz(2)})
    %keep a record of ind, which speeds up similar sizes
    [x,y]=meshgrid(1:size(data,2),1:size(data,1));
    z=mod(x+y-2,size(data,1))+1;
    ind=sub2ind(size(data),z,x);
    if prod(sz)<=1000 %to prevent a memory-hog, only keep ind for small sizes
        LUT{sz(1),sz(2)}=ind;
    end
else
    ind=LUT{sz(1),sz(2)};
end
data=data(ind);
end
function data=ComputeNonCryptHash_add_salt(data,UseVersion)
%Apply a few transformations to the hash to increase the spread.
%If this function is not added, the hashes of -12345678 and 12345678 will be very similar.
%A true salt would be to append the salt to the data, so this is not actually a salt.
saltHashLength=16*numel(data);
%Avoid an infinite recursion by using a fourth input:
salt=ComputeNonCryptHash(data,saltHashLength,'-v1',[]);
salt=ComputeNonCryptHash_shuffle_uint16_inv(salt);
if UseVersion>1
    salt=salt(end:-1:1);
end
data=mod(double(data).*double(salt),1+2^16);
data=uint16(data);
end
function hash=ComputeNonCryptHash_dec2hex(data)
%Look up the precomputed dec2hex for faster conversion.
persistent LUT
if isempty(LUT)
    LUT=upper(dec2hex(0:(-1+2^16),4));%Even though the default should already upper case.
end
data=double(data)+1;% plus() is not defined for uint16 on ML6.5
hash=LUT(data,:);
end
function data=ComputeNonCryptHash_logical_to_uint16(data)
if mod(numel(data),16) %Pad to 16 bits with 0.
    data(16*ceil(numel(data)/16))=0;
end
vector=uint16(2.^(15:-1:0))';
data=uint16(reshape(data,16,[]));
try
    data=data.*vector;
catch %No implicit expansion.
    %times() is not defined for uint16 on ML6.5
    data=double(data).*repmat(double(vector),[1 size(data,2)]);
    data=uint16(data);
end
data=uint16(sum(data,1)).';
end
function data=ComputeNonCryptHash_shuffle_uint16(data)
%Input should be uint16.
base=65537;%base=(1+(2^16));
key=479001600;%key=1*2*3*4*5*6*7*8*9*10*11*12;
data = uint16(mod(double(data) * key , base));
end
function data=ComputeNonCryptHash_shuffle_uint16_inv(data)
base=65537;%base=(1+(2^16));
%key=1*2*3*4*5*6*7*8*9*10*11*12;
% %Solution suggested by John D'Errico, https://www.mathworks.com/matlabcentral/answers/81859
% [G,C]=gcd(key,base);invKey=mod(C,base);
invKey=1919;
data=uint16(mod(double(data) * invKey,base));
end
function data=ComputeNonCryptHash_uint16_to_logical(data)
%uint16 Nx1 vector in, logical 16xN array out
persistent LUT
if isempty(LUT)
    LUT=dec2bin(0:(-1+2^16))=='1';
    LUT=LUT.';
end
data=double(data)+1;% plus() is not defined for uint16 on ML6.5
data=LUT(:,data);
end
function varargout=debug_hook(varargin)
%This function can be used to return several outputs (including warnings and errors), determined by
%the global variable HJW___test_suit___debug_hook_data.
%
%Every iteration the first element is removed from HJW___test_suit___debug_hook_data.
%
%When HJW___test_suit___debug_hook_data is empty or when returning a warning this functions returns
%the input unchanged.
global HJW___test_suit___debug_hook_data
if isempty(HJW___test_suit___debug_hook_data)
    varargout=varargin;return
end

element=HJW___test_suit___debug_hook_data(1);
HJW___test_suit___debug_hook_data(1)=[];

switch element.action
    case 'return'
        varargout=element.data;
    case 'warning'
        varargout=varargin;
        warning(element.data{:})
    case 'error'
        error(element.data{:})
    case 'warning_'
        varargout=varargin;
        warning_(element.data{:})
    case 'error_'
        error_(element.data{:})
end
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
function error_(options,varargin)
%Print an error to the command window, a file and/or the String property of an object.
%The error will first be written to the file and object before being actually thrown.
%
%The intention is to allow replacement of every error(___) call with error_(options,___).
%
% NB: the error trace that is written to a file or object may differ from the trace displayed by
% calling the builtin error function. This was only observed when evaluating code sections. 
%
%options.fid.boolean: if true print error to file (options.fid.fid)
%options.obj.boolean: if true print error to object (options.obj.obj)
%
%syntax:
%  error_(options,msg)
%  error_(options,msg,A1,...,An)
%  error_(options,id,msg)
%  error_(options,id,msg,A1,...,An)
%  error_(options,ME)               %equivalent to rethrow(ME)

%Parse input to find id, msg, stack and the trace str.
if isempty(options),options=struct;end%allow empty input to revert to default
if ~isfield(options,'fid'),options.fid.boolean=false;end
if ~isfield(options,'obj'),options.obj.boolean=false;end
if nargin==2
    %  error_(options,msg)
    %  error_(options,ME)
    if isa(varargin{1},'struct') || isa(varargin{1},'MException')
        ME=varargin{1};
        try
            stack=ME.stack;%Use the original call stack if possible.
            trace=get_trace(0,stack);
        catch
            [trace,stack]=get_trace(2);
        end
        id=ME.identifier;
        msg=ME.message;
        pat='Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback(';
        %This pattern may occur when using try error(id,msg),catch,ME=lasterror;end instead of
        %catching the MException with try error(id,msg),catch ME,end.
        %This behavior is not stable enough to robustly check for it, but it only occurs with
        %lasterror, so we can use that.
        if isa(ME,'struct') && numel(msg)>numel(pat) && strcmp(pat,msg(1:numel(pat)))
            %Strip the first line (which states 'error in function (line)', instead of only msg).
            msg(1:find(msg==10,1))='';
        end
    else
        [trace,stack]=get_trace(2);
        [id,msg]=deal('',varargin{1});
    end
else
    [trace,stack]=get_trace(2);
    if ~isempty(strfind(varargin{1},'%')) %The id can't contain a percent symbol.
        %  error_(options,msg,A1,...,An)
        id='';
        A1_An=varargin(2:end);
        msg=sprintf(varargin{1},A1_An{:});
    else
        %  error_(options,id,msg)
        %  error_(options,id,msg,A1,...,An)
        id=varargin{1};
        msg=varargin{2};
        if nargin>3
            A1_An=varargin(3:end);
            msg=sprintf(msg,A1_An{:});
        end
    end
end
ME=struct('identifier',id,'message',msg,'stack',stack);

%Print to object.
if options.obj.boolean
    msg_=msg;while msg_(end)==10,msg_(end)='';end%Crop trailing newline.
    if any(msg_==10)  % Parse to cellstr and prepend error.
        msg_=regexp_outkeys(['Error: ' msg_],char(10),'split'); %#ok<CHARTEN>
    else              % Only prepend error.
        msg_=['Error: ' msg_];
    end
    set(options.obj.obj,'String',msg_)
end

%Print to file.
if options.fid.boolean
    fprintf(options.fid.fid,'Error: %s\n%s',msg,trace);
end

%Actually throw the error.
rethrow(ME)
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
function [str,stack]=get_trace(skip_layers,stack)
if nargin==0,skip_layers=1;end
if nargin<2, stack=dbstack;end
stack(1:skip_layers)=[];

%Parse the ML6.5 style of dbstack (the name field includes full file location).
if ~isfield(stack,'file')
    for n=1:numel(stack)
        tmp=stack(n).name;
        if strcmp(tmp(end),')')
            %Internal function.
            ind=strfind(tmp,'(');
            name=tmp( (ind(end)+1):(end-1) );
            file=tmp(1:(ind(end)-2));
        else
            file=tmp;
            [ignore,name]=fileparts(tmp); %#ok<ASGLU>
        end
        [ignore,stack(n).file]=fileparts(file); %#ok<ASGLU>
        stack(n).name=name;
    end
end

%Parse Octave style of dbstack (the file field includes full file location).
persistent IsOctave,if isempty(IsOctave),IsOctave=exist('OCTAVE_VERSION', 'builtin');end
if IsOctave
    for n=1:numel(stack)
        [ignore,stack(n).file]=fileparts(stack(n).file); %#ok<ASGLU>
    end
end

%Create the char array with a (potentially) modified stack.
s=stack;
c1='>';
str=cell(1,numel(s)-1);
for n=1:numel(s)
    [ignore_path,s(n).file,ignore_ext]=fileparts(s(n).file); %#ok<ASGLU>
    if n==numel(s),s(n).file='';end
    if strcmp(s(n).file,s(n).name),s(n).file='';end
    if ~isempty(s(n).file),s(n).file=[s(n).file '>'];end
    str{n}=sprintf('%c In %s%s (line %d)\n',c1,s(n).file,s(n).name,s(n).line);
    c1=' ';
end
str=horzcat(str{:});
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
if debug_hook(false),EntryFunctionNames=ListLocalFunctions(str);end % Overwrite EntryFunctionNames.
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
function out=PatternReplace(in,pattern,rep)
%Functionally equivalent to strrep, but extended to more data types.
out=in(:)';
if numel(pattern)==0
    L=false(size(in));
elseif numel(rep)>numel(pattern)
    error('not implemented (padding required)')
else
    L=true(size(in));
    for n=1:numel(pattern)
        k=find(in==pattern(n));
        k=k-n+1;k(k<1)=[];
        %Now k contains the indices of the beginning of each match.
        L2=false(size(L));L2(k)=true;
        L= L & L2;
        if ~any(L),break,end
    end
end
k=find(L);
if ~isempty(k)
    for n=1:numel(rep)
        out(k+n-1)=rep(n);
    end
    if numel(rep)==0,n=0;end
    if numel(pattern)>n
        k=k(:);%Enforce direction.
        remove=(n+1):numel(pattern);
        idx=bsxfun_plus(k,remove-1);
        out(idx(:))=[];
    end
end
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
function [bit64,last16bits]=typecast_double_uint64(FP)
%Turn a double into a uint64 with the same binary representation.
%This is similar to typecast(FP,'uint64'); the difference being that this function ignores
%endianness, supports array inputs, and is slower.
%
%Because of missing support for some operations in ML6.5, this function returns the uint64 as a
%double. Because this may cause rounding errors, the last 16 bits are returned separately.

[M,E]=log2(FP);
signBit =-floor(sign(FP)/2-0.5);
exponent=E+1022;
mantissa=abs(M)*2-1;

%no plus() for integer types in ML6.5, so we need to use double, instead of uint64
bit64=zeros(size(FP));
bit64=bit64+(signBit*2^63);
bit64=bit64+(exponent*2^52);
bit64=bit64+(mantissa*2^52);
last16bits=mod(mantissa*2^52,2^16);

%correct the 0 and hard-code the special cases
L=isinf(FP);
bit64(FP==0)=0;
bit64(isnan(FP))=18444492273895866368;
bit64(L & FP>0)=9218868437227405312;%positive inf
bit64(L & FP<0)=18442240474082181120;%negative inf
last16bits(FP==0)=0;
last16bits(isnan(FP))=0;
last16bits(L)=0;
end
function str=unicode_to_char(unicode,encode_as_UTF16)
%Encode Unicode code points with UTF-16 on Matlab and UTF-8 on Octave.
%Input is implicitly converted to a row-vector.

persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end
if nargin==1
    encode_as_UTF16=~isOctave;
end
if encode_as_UTF16
    if all(unicode<65536)
        str=uint16(unicode);
        str=reshape(str,1,numel(str));%Convert explicitly to a row-vector.
    else
        %Encode as UTF-16.
        [char_list,ignore,positions]=unique(unicode); %#ok<ASGLU>
        str=cell(1,numel(unicode));
        for n=1:numel(char_list)
            str_element=unicode_to_UTF16(char_list(n));
            str_element=uint16(str_element);
            str(positions==n)={str_element};
        end
        str=cell2mat(str);
    end
    if ~isOctave
        str=char(str);%Conversion to char could trigger a conversion range error in Octave.
    end
else
    if all(unicode<128)
        str=char(unicode);
        str=reshape(str,1,numel(str));%Convert explicitly to a row-vector.
    else
        %Encode as UTF-8
        [char_list,ignore,positions]=unique(unicode); %#ok<ASGLU>
        str=cell(1,numel(unicode));
        for n=1:numel(char_list)
            str_element=unicode_to_UTF8(char_list(n));
            str_element=uint8(str_element);
            str(positions==n)={str_element};
        end
        str=cell2mat(str);
        str=char(str);
    end
end
end
function str=unicode_to_UTF16(unicode)
%Convert a single character to UTF-16 bytes.
%
%The value of the input is converted to binary and padded with 0 bits at the front of the string to
%fill all 'x' positions in the scheme.
%See https://en.wikipedia.org/wiki/UTF-16
%
% 1 word (U+0000 to U+D7FF and U+E000 to U+FFFF):
%  xxxxxxxx_xxxxxxxx
% 2 words (U+10000 to U+10FFFF):
%  110110xx_xxxxxxxx 110111xx_xxxxxxxx
if unicode<65536
    str=unicode;return
end
U=double(unicode)-65536;%Convert to double for ML6.5.
U=dec2bin(U,20);
str=bin2dec(['110110' U(1:10);'110111' U(11:20)]).';
end
function str=unicode_to_UTF8(unicode)
%Convert a single character to UTF-8 bytes.
%
%The value of the input is converted to binary and padded with 0 bits at the front of the string to
%fill all 'x' positions in the scheme.
%See https://en.wikipedia.org/wiki/UTF-8
if unicode<128
    str=unicode;return
end
persistent pers
if isempty(pers)
    pers=struct;
    pers.limits.lower=hex2dec({'0000','0080','0800', '10000'});
    pers.limits.upper=hex2dec({'007F','07FF','FFFF','10FFFF'});
    pers.scheme{2}='110xxxxx10xxxxxx';
    pers.scheme{2}=reshape(pers.scheme{2}.',8,2);
    pers.scheme{3}='1110xxxx10xxxxxx10xxxxxx';
    pers.scheme{3}=reshape(pers.scheme{3}.',8,3);
    pers.scheme{4}='11110xxx10xxxxxx10xxxxxx10xxxxxx';
    pers.scheme{4}=reshape(pers.scheme{4}.',8,4);
    for b=2:4
        pers.scheme_pos{b}=find(pers.scheme{b}=='x');
        pers.bits(b)=numel(pers.scheme_pos{b});
    end
end
bytes=find(pers.limits.lower<unicode & unicode<pers.limits.upper);
str=pers.scheme{bytes};
scheme_pos=pers.scheme_pos{bytes};
b=dec2bin(unicode,pers.bits(bytes));
str(scheme_pos)=b;
str=bin2dec(str.').';
end
function [unicode,isUTF8,assumed_UTF8]=UTF8_to_unicode(UTF8,print_to)
%Convert UTF-8 to the code points stored as uint32
%Plane 16 goes up to 10FFFF, so anything larger than uint16 will be able to hold every code point.
%
%If there a second output argument, this function will not return an error if there are encoding
%error. The second output will contain the attempted conversion, while the first output will
%contain the original input converted to uint32.
%
%The second input can be used to also print the error to a GUI element or to a text file.
if nargin<2,print_to=[];end
return_on_error= nargout==1 ;

UTF8=uint32(UTF8);
[assumed_UTF8,flag,ME]=UTF8_to_unicode_internal(UTF8,return_on_error);
if strcmp(flag,'success')
    isUTF8=true;
    unicode=assumed_UTF8;
elseif strcmp(flag,'error')
    isUTF8=false;
    if return_on_error
        error_(print_to,ME)
    end
    unicode=UTF8;%Return input unchanged (apart from casting to uint32).
end
end
function [UTF8,flag,ME]=UTF8_to_unicode_internal(UTF8,return_on_error)

flag='success';
ME=struct('identifier','HJW:UTF8_to_unicode:notUTF8','message','Input is not UTF-8.');

persistent isOctave,if isempty(isOctave),isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;end

if any(UTF8>255)
    flag='error';
    if return_on_error,return,end
elseif all(UTF8<128)
    return
end

for bytes=4:-1:2
    val=bin2dec([repmat('1',1,bytes) repmat('0',1,8-bytes)]);
    multibyte=UTF8>=val & UTF8<256;%Exclude the already converted chars.
    if any(multibyte)
        multibyte=find(multibyte);multibyte=multibyte(:).';
        if numel(UTF8)<(max(multibyte)+bytes-1)
            flag='error';
            if return_on_error,return,end
            multibyte( (multibyte+bytes-1)>numel(UTF8) )=[];
        end
        if ~isempty(multibyte)
            idx=bsxfun_plus(multibyte , (0:(bytes-1)).' );
            idx=idx.';
            multibyte=UTF8(idx);
        end
    else
        multibyte=[];
    end
    header_bits=[repmat('1',1,bytes-1) repmat('10',1,bytes)];
    header_locs=unique([1:(bytes+1) 1:8:(8*bytes) 2:8:(8*bytes)]);
    if numel(multibyte)>0
        multibyte=unique(multibyte,'rows');
        S2=mat2cell(multibyte,ones(size(multibyte,1),1),bytes);
        for n=1:numel(S2)
            bin=dec2bin(double(S2{n}))';
            %To view the binary data, you can use this: bin=bin(:)';
            %Remove binary header (3 byte example):
            %1110xxxx10xxxxxx10xxxxxx
            %    xxxx  xxxxxx  xxxxxx
            if ~strcmp(header_bits,bin(header_locs))
                %Check if the byte headers match the UTF-8 standard.
                flag='error';
                if return_on_error,return,end
                continue %leave unencoded
            end
            bin(header_locs)='';
            if ~isOctave
                S3=uint32(bin2dec(bin  ));
            else
                S3=uint32(bin2dec(bin.'));%Octave needs an extra transpose.
            end
            %Perform actual replacement.
            UTF8=PatternReplace(UTF8,S2{n},S3);
        end
    end
end
end
function warning_(options,varargin)
%Print a warning to the command window, a file and/or the String property of an object.
%The lastwarn state will be set if the warning isn't thrown with warning().
%The printed call trace omits this function, but the warning() call does not.
%
%The intention is to allow replacement of most warning(___) call with warning_(options,___). This
%does not apply to calls that query or set the warning state.
%
%options.con:         if true print warning to command window with warning()
%options.fid.boolean: if true print warning to file (options.fid.fid)
%options.obj.boolean: if true print warning to object (options.obj.obj)
%
%syntax:
%  warning_(options,msg)
%  warning_(options,msg,A1,...,An)
%  warning_(options,id,msg)
%  warning_(options,id,msg,A1,...,An)

if isempty(options),options=struct;end%Allow empty input to revert to default.
if ~isfield(options,'con'),options.con=false;end
if ~isfield(options,'fid'),options.fid.boolean=false;end
if ~isfield(options,'obj'),options.obj.boolean=false;end
if nargin==2 || ~isempty(strfind(varargin{1},'%'))%The id can't contain a percent symbol.
    %  warning_(options,msg,A1,...,An)
    [id,msg]=deal('',varargin{1});
    if nargin>3
        A1_An=varargin(2:end);
        msg=sprintf(msg,A1_An{:});
    end
else
    %  warning_(options,id,msg)
    %  warning_(options,id,msg,A1,...,An)
    [id,msg]=deal(varargin{1},varargin{2});
    if nargin>3
        A1_An=varargin(3:end);
        msg=sprintf(msg,A1_An{:});
    end
end

if options.con
    if ~isempty(id)
        warning(id,'%s',msg)
    else
        warning(msg)
    end
else
    if ~isempty(id)
        lastwarn(msg,id);
    else
        lastwarn(msg)
    end
end

if options.obj.boolean
    msg_=msg;while msg_(end)==10,msg_(end)=[];end%Crop trailing newline.
    if any(msg_==10)  % Parse to cellstr and prepend warning.
        msg_=regexp_outkeys(['Warning: ' msg_],char(10),'split'); %#ok<CHARTEN>
    else              % Only prepend warning.
        msg_=['Warning: ' msg_];
    end
    set(options.obj.obj,'String',msg_)
end

if options.fid.boolean
    skip_layers=2;%Remove this function and the get_trace function from the trace.
    trace=get_trace(skip_layers);
    fprintf(options.fid.fid,'Warning: %s\n%s',msg,trace);
end
end