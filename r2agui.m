function varargout = r2agui(varargin)
% R2AGUI Application M-file for r2agui.fig
%
%   This program converts the Philips PAR and REC image file format
%   to the Analu=yze format used by SPM and other image processing packages
%   CAUTION: works for the very recent PAR & REC file version V4 only! For
%   the older V3 PAR file format conversion see
%   http://www.fss.uu.nl/psn/pionier/downloads.
%   Type 'r2agui' on the matlab command line to invoke the graphical user
%   interface. Make sure that r2agui.m and r2agui.fig are both in the
%   matlab search path.
%
%   Created by:
%   Erno Hermans and Bas Neggers, Helmholtz Institute,
%   Utrecht University/Dept of Brain Research, University Medical Center
%
%   Version 2.1, released 31-3-2006

%
% Changed for v4.1 & v4.2 using read_par.m from r2agui 2.4 : Jed Dobson
% 07/15/2008
%
% Modified by Yu-CHien Wu 2009/07/21
% in read_par file add tech (MR technique), diffusion (diffusion flag),
% bvalno (b value number) and bvecno (diffusion gradient number)
%

% added by Jed
warning off MATLAB:divideByZero

if nargin == 0  % LAUNCH GUI

    fig = openfig(mfilename,'reuse');

    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    guidata(fig, handles);

    if nargout > 0
        varargout{1} = fig;
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr);
    end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and
%| sets objects' callback properties to call them through the FEVAL
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = batchtoggle_Callback(h, eventdata, handles, varargin)

switch(get(handles.batchtoggle,'value')),
    case 1
        disp('Process all volumes in directory');
        set(handles.batchtoggle,'string','Batch on');
    case 0
        disp('Process selected volume only');
        set(handles.batchtoggle,'string','Batch off');
end

% --------------------------------------------------------------------
function varargout = openpar_Callback(h, eventdata, handles, varargin)
global pathpar
global file
if pathpar==0,
    pathpar=[];
end
[file, pathpar] = uigetfile([pathpar,'*.par;',pathpar,'*.PAR'], 'Select PAR file');

parfile = [pathpar, file];
Pars=read_par(parfile);

set(handles.voxelsizex, 'String',Pars.vox(1));
set(handles.voxelsizey, 'String',Pars.vox(2));
set(handles.voxelsizez, 'String',Pars.vox(3));
set(handles.dynamics, 'String',Pars.dyn);
set(handles.RT_Version, 'String',Pars.ResToolsVersion);
set(handles.filenaam, 'String',parfile);
prefix = get(handles.prefix, 'String');

subaan = get(handles.subdir, 'Value');
if subaan ==1

    subdir = [prefix,file(1:(length(file)-4))];
else
    subdir = '';
end
set(handles.outputfile, 'String',fullfile(pathpar,subdir, [prefix,file(1:(length(file)-4)),'-001.hdr']));

set(handles.convert, 'Enable', 'on')


% --------------------------------------------------------------------
function varargout = voxelsizex_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = voxelsizey_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = voxelsizez_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = dynamics_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = filenaam_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = prefix_Callback(h, eventdata, handles, varargin)
setfoldernames(handles);

% --------------------------------------------------------------------
function varargout = subdir_Callback(h, eventdata, handles, varargin)
setfoldernames(handles);

% --------------------------------------------------------------------
function varargout = usealtfolder_Callback(h, eventdata, handles, varargin)
setfoldernames(handles);

% --------------------------------------------------------------------
function varargout = altfolder_Callback(h, eventdata, handles, varargin)
setfoldernames(handles);

function varargout=setfoldernames(handles)
global pathpar
global file
if ispc,
    WinLinSlash='\';
else
    WinLinSlash='/';
end
prefix = get(handles.prefix, 'String');
subaan = get(handles.subdir, 'Value');
altfolderaan=get(handles.usealtfolder, 'Value');

if subaan ==1;

    subdir = [prefix,file(1:(length(file)-4))];
else
    subdir = '';
end

if altfolderaan ==1;
    outfoldername = get(handles.altfolder, 'String');
    if ~strcmp(outfoldername(end),WinLinSlash) %check if foldername has a trailing slash, add when necessary
        outfoldername=[outfoldername,WinLinSlash];
        set(handles.altfolder, 'String',outfoldername);
    end
else
    outfoldername = pathpar;
end


set(handles.outputfile, 'String',fullfile(outfoldername,subdir, [prefix,file(1:(length(file)-4)),'-001.hdr']));


% --------------------------------------------------------------------
function varargout = edit8_Callback(h, eventdata, handles, varargin)







% --------------------------------------------------------------------r2a

function varargout = convert_Callback(h, eventdata, handles, varargin)
global file;
global pathpar;
if ispc,
    WinLinSlash='\';
else
    WinLinSlash='/';
end
if get(handles.batchtoggle,'value')==0
    no_files = 1;
    filelist{1}=file;
else
    if ispc,    %check whether this is windows
        dl=dir([pathpar,'*.par']);
        no_files = length(dl);
        for i=1:length(dl),
            filelist{i}=dl(i).name;
        end
    else
        dlL=dir([pathpar,'*.par']);
        dlU=dir([pathpar,'*.PAR']);
        no_files = length(dlL)+length(dlU);
        for i=1:length(dlL),
            filelist{i}=dlL(i).name;
        end
        if length(dlL)==0,
            i=0;
        end
        for i2=1:length(dlU),
            filelist{i+i2}=dlU(i2).name;
        end
    end
end
for i=1 : no_files(1),
    parfile = [pathpar,filelist{i}];


    prefix = get(handles.prefix, 'String');
    Parameters=read_par(parfile);
    if Parameters.problemreading==1,
        disp(['Skipping volume ',parfile,' because of reading errors.']);
    else
        subaan = get(handles.subdir, 'Value');
        altfolderaan=get(handles.usealtfolder, 'Value');
        if subaan ==1

            subdir = [prefix,filelist{i}(1:(length(filelist{i})-4))];
        else
            subdir = '';
        end

        if altfolderaan ==1;
            outfoldername = get(handles.altfolder, 'String');
        else
            outfoldername = pathpar;
        end


        set(handles.outputfile, 'String',fullfile(outfoldername,subdir, [prefix,filelist{i}(1:(length(file)-4)),'-001.hdr']));

        set(handles.convert, 'Enable', 'on')

        Vox = Parameters.vox;

        Recfile=filelist{i};
        if strcmp(Recfile(end-2:end),'par')
            Recfile(end-2:end)='rec';
        elseif strcmp(Recfile(end-2:end),'PAR')
            Recfile(end-2:end)='REC';
        end
        if subaan ==1,

            subdir = [prefix,filelist{i}(1:(length(filelist{i})-4))];
            if exist([outfoldername,subdir])<7
                mkdir(outfoldername,subdir);
            end


        else
            subdir = '';
        end


        absDir = [outfoldername subdir];
        %   mkdir (DirName, filename );

        % converting  par-file parameters into  analyze header parameters

        Precision = strcat ('int', Parameters.bit);
        Precision = char (Precision);
        Size = Parameters.dim(1)*Parameters.dim(2)*Parameters.dim(3);
        SizeSlice = Parameters.dim(1)*Parameters.dim(2);
        ID1 = fopen ([pathpar,Recfile],'r','l');
        %Vox = Parameter.fov ./ Parameter.dim;
        %Vox = input('Type Voxel size (mm):');
        %if length(Vox)~=3,
        %    disp('Wrong voxel size input');
        %    break;
        %end

        Dim = Parameters.dim;
        Type = spm_type (Precision);
        Offset = 0;
        Descrip = char (Parameters.name);
        Scale = 1;
        Orign = Parameters.fov ./Vox /2 ;
        BytesPerValue=str2num(Parameters.bit)/8;



        switch(Parameters.ResToolsVersion)
            case 'V3'
                % old method of reading & writing files for V3 data
                disp ([' Start to convert scan:' Recfile ]);
                for j =1 : size(Paramaters.slice_index,1)

                    % default ordering of slices/dynamics changed between V3 and
                    % V4. NOTE: in near future slice info will be read from PAR
                    % file, and all versions should work right away.
                    switch(Parameters.slicessorted)
                        case 1 % in version 3, data was by default ordered dyn 1, slice 1, dyn 1, slice 2, etc
                            Data = fread (ID1,Size, Precision);
                            Inputvolume  = zeros(Parameters.dim(2),Parameters.dim(1),Parameters.dim(3));
                            Inputvolume = reshape(Data,Parameters.dim(2),Parameters.dim(1),Parameters.dim(3));
                        case 2 % in version 4, data is by default ordered dyn 1, slice 1, dyn 2, slice 1, etc
                            clear InputVolume;
                            Inputvolume  = zeros(Parameters.dim(2),Parameters.dim(1),Parameters.dim(3));
                            for slice=1:Parameters.dim(3),
                                fseek(ID1,(j-1+Parameters.dyn*(slice-1))*SizeSlice*BytesPerValue,-1);
                                InputVolume(:,:,slice) = fread (ID1,SizeSlice, Precision);
                            end
                            Inputvolume = reshape(InputVolume,Parameters.dim(2),Parameters.dim(1),Parameters.dim(3));
                    end
                    VolName =  [absDir WinLinSlash subdir '-' num2str(j,'%03.0f') '.img'];

                    ID2 = fopen (VolName, 'w');
                    fwrite (ID2, Inputvolume (:,:,:),Precision);
                    fclose (ID2);

                    P = VolName;
                    spm_hwrite (P,Dim,Vox,Scale,Type ,Offset,round(Orign),Descrip);

                    disp (['Write file: ' subdir '-' num2str(j,'%03.0f')]);
                end
            case 'V4'
                % new: loop slices (as in slice_index) and open and close files
                % along the way (according to info in index on dynamic and mr_type)
                iSlice=Parameters.slice_index;
                iSlice(:,12)=iSlice(:,8).*iSlice(:,10).*iSlice(:,11)/8; % add column containing size of slice in bytes
                order_slices=iSlice(:,7); % get order in which slices are stored in REC file.
                [os,i]=sort(order_slices); % sort them
                bytespslice_sorted=iSlice(i,12); %sort bytes per slice info acc. to index in REC file
                fileposSlice_sorted=[cumsum(bytespslice_sorted)]; % sum bytes per slice cumulatively
                fileposSlice_sorted=[0;fileposSlice_sorted(1:end-1)]; % start in file in bytes of each slice (sorted)
                index_orig=[1:size(order_slices)];
                fileposSlice=fileposSlice_sorted(index_orig(i)); % unsort file position of slice in bytes.
                iSlice(:,13)=fileposSlice; % add column containing start position in bytes of this slice in the file
                
                %now sort entire slice_index according to dynamics
                %(fastest varying) and mr_type parameters (slowest varying)
                iSlices_sorted = sortrows(iSlice,[6 3 1]);
                
                nLine=0;
                NewFile=[1; (diff(iSlices_sorted(:,3))~=0 | diff(iSlices_sorted(:,6))~=0)]; % determine whether new file has to be opened (change of dynamic/mr_type)
                nr_mrtypes=length(unique(iSlices_sorted(:,6))); % determine number of interleaved image types (e.g. angio)
                while nLine<size(iSlices_sorted,1);
                    
                    nLine=nLine+1;
                    if NewFile(nLine)
                       %close previous file
                       if nLine>1
                           fclose(ID2);
                       end
                       if nr_mrtypes>1
                           mrtype_suffix=num2str(iSlices_sorted(nLine,6),'-%03.0i');
                       else
                           mrtype_suffix='';
                       end
                       Precision = strcat ('int', num2str(iSlices_sorted(nLine,8)));
                       Precision = char (Precision);
                       cDim=[iSlices_sorted(nLine,10:11),Dim(3)];
                       VolName =  [absDir WinLinSlash subdir mrtype_suffix '-' num2str(iSlices_sorted(nLine,3),'%03i') '.img'];
                       ID2 = fopen (VolName, 'w');
                       P = VolName;
                       spm_hwrite (P,cDim,Vox,Scale,Type ,Offset,round(Orign),Descrip);
                       disp (['Write file: ' subdir '-' num2str(iSlices_sorted(nLine,3),'%03.0f')]);                         
                    end
                       
                   fseek(ID1,iSlices_sorted(nLine,13),-1);
                   SliceData = fread (ID1,cDim(1)*cDim(2), Precision);
                   fwrite(ID2, SliceData,Precision);
                    
                end
                fclose(ID2);
                fclose(ID1);
                
                
                
            otherwise
                disp(['Sorry, but data format extracted using Philips Research File format ', ...
                    Parameters.ResToolsVersion,' was not known at the time the r2agui software was developed']);
        end
    end
end


function  [WorkDir, filename]  = l_split(P,ext)

n = size (P);
n = n(1);
if n == 0
    arrayout = 'error'
end

if n > 0
    filename = P(1,:);

    ind = findstr(filename,WinLinSlash);

    if ~isempty(ind)

        WorkDir = filename(1: ind (length (ind))); %extract WorkDir

        ind2 = findstr( filename, ext);
        if ~isempty(ind2)
            filename = filename (ind (length (ind))+1 :ind2-1);
        end

    end
end



function  [par]  = read_par (parfile);

%read version number of Philips research tools
%Research tools are used to extract data from database; dataformats differ considerably
%between versions. R2AGUI now handles V3 and V4

par.problemreading=0; % will be set to 1 when an error occurs
parameter = textread (parfile,'%s',8, 'headerlines',7);
par.ResToolsVersion=parameter{8};

switch(par.ResToolsVersion)
    case 'V3'
        parameter = textread (parfile,'%s',5, 'delimiter', '.:','headerlines',11);
        par.name = parameter (4);
        % jed add
        test = char(par.name);
        if strcmp(test,'Sense Ref')
           par.dyn = 1;
           par.slice = 1;
           par.vox = 1;
           return;
        end
        parameter = textread (parfile,'%u', 16,'delimiter','.:Acquisitionr','headerlines',15);
        par.scno = parameter (16);
        parameter = textread (parfile,'%u', 31,'delimiter','.:Max.numberofslices/locations','headerlines',20);
        par.slice = parameter (31);
        parameter = textread (parfile,'%u', 23,'delimiter','.:Max.numberofdynamics','headerlines',21);
        par.dyn = parameter (23);
        parameter = textread (parfile,'%s', 25,'delimiter','.:Imagepixelsize[orbits]','headerlines',23);
        par.bit = parameter {25};
        parameter = textread (parfile,'%u', 24,'delimiter','.:Reconresolution(x,y)','headerlines',28);
        %parameter = textread (parfile,'%u', 24,'delimiter','.:Scanresolution(x,y)','headerlines',26);
        x = parameter (23);
        y = parameter (24);
        z = par.slice;
        par.dim = [x,y,z];
        parameter = textread (parfile,'%s', 20,'headerlines',90);
        par.sliceorient  = parameter (20);
        parameter = textread (parfile,'%f', 22,'delimiter','.:FOV(ap,fh,rl)[mm]','headerlines',31);

        if strcmp (par.sliceorient, '1')  ; % slice orientation: transversal
            fovx = parameter (20);
            fovy = parameter (22);
            par.sliceorient;
        end;

        if strcmp (par.sliceorient, '2')   ;% slice orientation: sagital
            fovx = parameter (21)
            fovy = parameter (20)
            par.sliceorient
        end;

        if strcmp (par.sliceorient,'3')   ;% slice orientation: coronal
            fovx = parameter (22)
            fovy = parameter (21)
            par.sliceorient
        end;

        parameter = textread (parfile,'%f', 21,'delimiter','.:Slicethickness[mm]','headerlines',32);
        par.slth = parameter (21);
        parameter = textread (parfile,'%f', 15,'delimiter','.:Slicegap[mm]','headerlines',33);
        par.gap= parameter (15);
        fovz = (par.gap + par.slth)*par.slice;
        par.fov = [fovx,fovy,fovz];
        parameter = textread (parfile,'%f', 39,'delimiter','.:Angulationmidslice(ap,fh,rl)[degr]','headerlines',35);
        par.angAP = parameter (37);
        par.angFH = parameter (38);
        par.angRL = parameter (39);
        parameter = textread (parfile,'%f', 36,'delimiter','.:OffCentremidslice(ap,fh,rl)[mm]','headerlines',36);
        par.offAP= parameter (34);
        par.offFH= parameter (35);
        par.offRL= parameter (36);
        parameter = textread (parfile,'%s',24, 'headerlines',88);
        voxx = str2num(parameter{23});
        voxy = str2num(parameter{24});
        voxz = par.slth + par.gap;
        par.vox=[voxx voxy voxz];
        parameternextline = textread (parfile,'%s',24, 'headerlines',90);
        if (parameternextline{1}-parameter{1})>0,
            par.slicessorted=1;
        else
            par.slicessorted=2;
        end

        clear parameter;

    case {'V4','V4.1','V4.2'}
        parameter = textread (parfile,'%s',5, 'delimiter', ':','headerlines',13);
        par.name = parameter (2); % for scantechnique in description/better of anonimity        
        % jed add
        test = char(par.name);
        if strcmp(test,'Sense Ref')
           par.dyn = 1;
           par.slice = 1;
           par.vox = 1;
           return
        end
	parameter = textread (parfile,'%u', 16,'delimiter','.:Acquisitionr','headerlines',16);
        par.scno = parameter (16);
        parameter = textread (parfile,'%u', 31,'delimiter','.:Max.numberofslices/locations','headerlines',21);
        par.slice = parameter (31);
        parameter = textread (parfile,'%u', 23,'delimiter','.:Max.numberofdynamics','headerlines',22);
        par.dyn = parameter (23);

        % read first six rows of indices from per-slice lines in PAR file
  	slice_index=textread (parfile,'','delimiter',' ','headerlines',90,'commentstyle','shell');

 	par.slice_index=slice_index(:,1:11);

        par.RT=(slice_index(end,32)-slice_index(1,32))/(par.dyn-1); % estimate scan-duration from dtime PAR file row
        %read first line of slice acquisition info        
	switch par.ResToolsVersion            
	case 'V4'                
		parameter = textread (parfile,'%s',41, 'headerlines',91); % V4
        case 'V4.1'                
		parameter = textread (parfile,'%s',41, 'headerlines',97); % PAR file format v4.1 is very similar to V4, only 5 extra diffusion columns are added to slice lines            
	case 'V4.2'                
	parameter = textread (parfile,'%s',41, 'headerlines',100); % PAR file format v4.2 is very similar to V4, only 5 extra diffusion columns are added to slice lines
        end
        if length(parameter)<15,
            disp(sprintf('Problem: No actual slices measured for volume %s.\nVolume might be corrupt.', parfile));
            par.problemreading=1;
        else
            par.sliceorient  = parameter{26};

            x = str2num(parameter{10});
            y = str2num(parameter{11});
            z = par.slice;
            par.dim = [x,y,z];

	    par.rescale_slope=str2num(parameter{13});
            par.rescale_interc=str2num(parameter{12});

            par.bit = parameter{8};
            par.slth = str2num(parameter{23});
            par.gap = str2num(parameter{24});
            voxx = str2num(parameter{29});
            voxy = str2num(parameter{30});
            voxz = par.slth+par.gap;
            par.vox=[voxx voxy voxz];

            % to check whether slices are ordered (ie if all slices are
            % written sequentially, or all volumes)
            parameternextline = textread (parfile,'%s',24, 'headerlines',92);
            if (parameternextline{1}-parameter{1})>0,
                par.slicessorted=1;
            else
                par.slicessorted=2;
            end

            parameter = textread (parfile,'%f', 22,'delimiter','.:FOV(ap,fh,rl)[mm]','headerlines',30);
            fovx = parameter (21);
            fovy = parameter (22);
            fovz = parameter (20);
            par.fov = [fovx,fovy,fovz]; 
            par.fov_apfhrl=[parameter(20) parameter(21) parameter(22)]; 
            parameter = textread (parfile,'%f', 39,'delimiter','.:Angulationmidslice(ap,fh,rl)[degr]','headerlines',32);
            fovx = parameter (21);
            fovy = parameter (22);
            fovz = parameter (20);
            par.fov = [fovx,fovy,fovz]; % leave this fov definition in only for backwards compatibility for V3 data
            par.fov_apfhrl=[parameter(20) parameter(21) parameter(22)]; %actually used now to calculate angulation etc

            parameter = textread (parfile,'%f', 39,'delimiter','.:Angulationmidslice(ap,fh,rl)[degr]','headerlines',32);
            par.angAP = parameter (37);
            par.angFH = parameter (38);
            par.angRL = parameter (39);
            parameter = textread (parfile,'%f', 36,'delimiter','.:OffCentremidslice(ap,fh,rl)[mm]','headerlines',33);

            par.offAP= parameter (34);
            par.offFH= parameter (35);
            par.offRL= parameter (36);

        end
        
        % added by Yu-Chien Wu 2009/07/21
        % read technique info:
        % works for V4.1, not sure about other versions.   Yu-Chien Wu   2009/07/21
        parameter = textread (parfile,'%s',2,'delimiter',':','headerlines',26);
        par.tech = parameter {2};
                
        % read diffusion parameters: 
        % works for V4.1, not sure about other versions.   Yu-Chien Wu   2009/07/21
        parameter = textread (parfile,'%s',10,'delimiter',':', 'headerlines',41);
        par.diffusion = parameter {2};
        par.bvalno = str2num(parameter{6});
        par.bvecno = str2num(parameter{8});        
end



% --------------------------------------------------------------------
function varargout = pushbutton3_Callback(h, eventdata, handles, varargin)



function [s] = spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP)
% writes a header
% (function copied from spm99, so spm99 does not have to be present)
% FORMAT [s] = spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);
%
% P       - filename 	     (e.g 'spm' or 'spm.img')
% DIM     - image size       [i j k [l]] (voxels)
% VOX     - voxel size       [x y z [t]] (mm [sec])
% SCALE   - scale factor
% TYPE    - datatype (integer - see spm_type)
% OFFSET  - offset (bytes)
% ORIGIN  - [i j k] of origin  (default = [0 0 0])
% DESCRIP - description string (default = 'spm compatible')
%
% s       - number of elements successfully written (should be 348)
%__________________________________________________________________________
%
% spm_hwrite writes variables from working memory into a SPM/ANALYZE
% compatible header file.  The 'originator' field of the ANALYZE format has
% been changed to ORIGIN in the SPM version of the header. funused1
% of the ANALYZE format is used for SCALE
%
% see also dbh.h (ANALYZE) spm_hread.m and spm_type.m
%
%__________________________________________________________________________
% @(#)spm_hwrite.m	2.2 99/10/29


% ensure correct suffix {.hdr} and open header file
%---------------------------------------------------------------------------
P               = P(P ~= ' ');
q    		= length(P);
if q>=4 & P(q - 3) == '.', P = P(1:(q - 4)); end;
P     		= [P '.hdr'];

% For byte swapped data-types, also swap the bytes around in the headers.
mach = 'native';
if spm_type(TYPE,'swapped'),
    if spm_platform('bigend'),
        mach = 'ieee-le';
    else,
        mach = 'ieee-be';
    end;
    TYPE = spm_type(spm_type(TYPE));
end;
fid             = fopen(P,'w',mach);

if (fid == -1),
    error(['Error opening ' P '. Check that you have write permission.']);
end;
%---------------------------------------------------------------------------
data_type 	= ['dsr      ' 0];

P     		= [P '                  '];
db_name		= [P(1:17) 0];

% set header variables
%---------------------------------------------------------------------------
DIM		= DIM(:)'; if size(DIM,2) < 4; DIM = [DIM 1]; end
VOX		= VOX(:)'; if size(VOX,2) < 4; VOX = [VOX 0]; end
dim		= [4 DIM(1:4) 0 0 0];
pixdim		= [0 VOX(1:4) 0 0 0];
vox_offset      = OFFSET;
funused1	= SCALE;
glmax		= 1;
glmin		= 0;
bitpix 		= 0;
descrip         = zeros(1,80);
aux_file        = ['none                   ' 0];
origin          = [0 0 0 0 0];

%---------------------------------------------------------------------------
if TYPE == 1;   bitpix = 1;  glmax = 1;        glmin = 0;	end
if TYPE == 2;   bitpix = 8;  glmax = 255;      glmin = 0;	end
if TYPE == 4;   bitpix = 16; glmax = 32767;    glmin = 0;  	end
if TYPE == 8;   bitpix = 32; glmax = (2^31-1); glmin = 0;	end
if TYPE == 16;  bitpix = 32; glmax = 1;        glmin = 0;	end
if TYPE == 64;  bitpix = 64; glmax = 1;        glmin = 0;	end

%---------------------------------------------------------------------------
if nargin >= 7; origin = [ORIGIN(:)' 0 0];  end
if nargin <  8; DESCRIP = 'spm compatible'; end

d          	= 1:min([length(DESCRIP) 79]);
descrip(d) 	= DESCRIP(d);

fseek(fid,0,'bof');

% write (struct) header_key
%---------------------------------------------------------------------------
fwrite(fid,348,		'int32');
fwrite(fid,data_type,	'char' );
fwrite(fid,db_name,	'char' );
fwrite(fid,0,		'int32');
fwrite(fid,0,		'int16');
fwrite(fid,'r',		'char' );
fwrite(fid,'0',		'char' );

% write (struct) image_dimension
%---------------------------------------------------------------------------
fseek(fid,40,'bof');

fwrite(fid,dim,		'int16');
fwrite(fid,'mm',	'char' );
fwrite(fid,0,		'char' );
fwrite(fid,0,		'char' );

fwrite(fid,zeros(1,8),	'char' );
fwrite(fid,0,		'int16');
fwrite(fid,TYPE,	'int16');
fwrite(fid,bitpix,	'int16');
fwrite(fid,0,		'int16');
fwrite(fid,pixdim,	'float');
fwrite(fid,vox_offset,	'float');
fwrite(fid,funused1,	'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'float');
fwrite(fid,0,		'int32');
fwrite(fid,0,		'int32');
fwrite(fid,glmax,	'int32');
fwrite(fid,glmin,	'int32');

% write (struct) image_dimension
%---------------------------------------------------------------------------
fwrite(fid,descrip,	'char');
fwrite(fid,aux_file,    'char');
fwrite(fid,0,           'char');
fwrite(fid,origin,      'int16');
if fwrite(fid,zeros(1,85), 'char')~=85
    fclose(fid);
    spm_unlink(P);
    error(['Error writing ' P '. Check your disk space.']);
end

s   = ftell(fid);
fclose(fid);

function T = spm_type(x, arg)
% translates data type specifiers between SPM & Matlab representations
% FORMAT T = spm_type(x, arg)
% x    - specifier
% T    - type
% arg  - optional string argument, can be
%	 - 'swapped' - if type is byteswapped return 1.
%	 - 'maxval'  - return maximum allowed value.
%	 - 'minval'  - return minimum allowed value.
%	 - 'nanrep'  - return 1 if there is a NaN representation.
%	 - 'bits'    - return the number of bits per voxel.
%	 - 'intt'    - return 1 if values rounded to nearest integer.
%_______________________________________________________________________
%
% Original format specifiers are based on ANALYZE.  If the input is
% a number then the corresponding matlab string is returned by default.
% If the input is a string then the appropriate TYPE is returned.
% However, if the optional arg argument is supplied then other
% information will be returned instead.
%
% With no arguments, a list of data types is returned.
%
% Additional support was added for signed bytes, unsigned short and
% unsigned int (by adding 128 to the format specifiers for unsigned bytes
% signed short and signed int).  Byte swapped datatypes have the same
% identifiers as the non-byte-swapped versions, multiplied by a factor of
% 256.
%_______________________________________________________________________
% @(#)spm_type.m	2.3 John Ashburner, Andrew Holmes 99/04/27


prec = str2mat('uint8','int16','int32','float','double','int8','uint16','uint32','uint8','int16','int32','float','double','int8','uint16','uint32');
types   = [    2      4      8   16   64   130    132    136,   512   1024   2048 4096 16384 33280  33792  34816];
swapped = [    0      0      0    0    0     0      0      0,     1      1      1    1     1     1      1      1];
maxval  = [2^8-1 2^15-1 2^31-1  Inf  Inf 2^7-1 2^16-1 2^32-1, 2^8-1 2^15-1 2^31-1  Inf   Inf 2^8-1 2^16-1 2^32-1];
minval  = [    0  -2^15  -2^31 -Inf -Inf  -2^7      0      0,     0  -2^15  -2^31 -Inf  -Inf  -2^7      0      0];
nanrep  = [    0      0      0    1    1     0      0      0,     0      0      0    1     1     0      0      0];
bits    = [    8     16     32   32   64     8     16     32,     8     16     32   32    64     8     16     32];
intt    = [    1      1      1    0    0     1      1      1,     1      1      1    0     0     1      1      1];

if nargin==0,
    T=types;
    return;
end;

if ischar(x),
    sel = [];
    msk = find(swapped==0);
    for i=msk,
        if strcmp(deblank(prec(i,:)),deblank(x)),
            sel = i;
            break;
        end;
    end;
else,
    sel = find(types == x);
end;
if nargin == 1,
    if ischar(x),
        if isempty(sel), T = NaN;
        else, T = types(sel); end;
    else,
        if isempty(sel), T = 'unknown';
        else, T = deblank(prec(sel,:)); end;
    end;
elseif isempty(sel),
    T = NaN;
else,
    switch lower(arg)
        case 'swapped', T = swapped(sel);
        case 'maxval',  T = maxval(sel);
        case 'minval',  T = minval(sel);
        case 'nanrep',  T = nanrep(sel);
        case 'bits',    T = bits(sel);
        case 'intt',    T = intt(sel);
        otherwise,      T = NaN;
    end;
end;


function varargout=spm_platform(varargin)
% Platform specific configuration parameters for SPM
%
% FORMAT ans = spm_platform(arg)
% arg  - optional string argument, can be
%        - 'bigend'  - return whether this architecture is bigendian
%                      - Inf - is not IEEE floating point
%                      - 0   - is little end
%                      - 1   - big end
%        - 'filesys' - type of filesystem
%                      - 'unx' - UNIX
%                      - 'win' - DOS
%                      - 'mac' - Macintosh
%                      - 'vms' - VMS
%        - 'sepchar' - returns directory separator
%        - 'rootlen' - returns number of chars in root directory name
%        - 'user'    - returns username
%        - 'tempdir' - returns name of temp directory
%
% FORMAT PlatFontNames = spm_platform('fonts')
% Returns structure with fields named after the generic (UNIX) fonts, the
% field containing the name of the platform specific font.
%
% FORMAT PlatFontName = spm_platform('font',GenFontName)
% Maps generic (UNIX) FontNames to platform specific FontNames
%
% FORMAT SPM_PLATFORM = spm_platform('init',comp)
% Initialises platform specific parameters in global SPM_PLATFORM
% (External gateway to init_platform(comp) subfunction)
% comp         - computer to use [defaults to MatLab's `computer`]
% SPM_PLATFORM - copy of global SPM_PLATFORM
%
% FORMAT spm_platform
% Initialises platform specific parameters in global SPM_PLATFORM
% (External gateway to init_platform(computer) subfunction)
%
% FORMAT spm_platform('clear')
% Clears global SPM_PLATFORM containing platform specific parameters
%
%                           ----------------
% SUBFUNCTIONS:
%
% FORMAT init_platform(comp)
% Initialise platform specific parameters in global SPM_PLATFORM
% comp         - computer to use [defaults to MatLab's `computer`]
%
%-----------------------------------------------------------------------
%
% Since calls to spm_platform will be made frequently, most platform
% specific parameters are stored as a structure in the global variable
% SPM_PLATFORM. Subsequent calls use the information from this global
% variable, if it exists.
%
% Platform specific difinitions are contained in the data structures at
% the beginning of the init_platform subfunction at the end of this
% file.
%_______________________________________________________________________
% @(#)spm_platform.m	2.10 Matthew Brett 00/11/08


%-Initialise
%-----------------------------------------------------------------------
global SPM_PLATFORM
if isempty(SPM_PLATFORM), init_platform, end

if nargin==0, return, end


switch lower(varargin{1}), case 'init'                  %-(re)initialise
    %=======================================================================
    init_platform(varargin{2:end})
    varargout = {SPM_PLATFORM};

    case 'clear'                                       %-Clear SPM_PLATFORM
        %=======================================================================
        clear global SPM_PLATFORM

    case 'bigend'                      %-Return endian for this architecture
        %=======================================================================
        varargout = {SPM_PLATFORM.bigend};
        if ~finite(SPM_PLATFORM.bigend),
            if isnan(SPM_PLATFORM.bigend)
                error(['I don''t know if "',computer,'" is big-endian.'])
            else
                error(['I don''t think that "',computer,...
                    '" uses IEEE floating point ops.'])
            end
        end

    case 'filesys'                                      %-Return file system
        %=======================================================================
        varargout = {SPM_PLATFORM.filesys};

    case 'sepchar'                         %-Return file separator character
        %=======================================================================
        warning('use filesep instead (supported by MathWorks)')
        varargout = {SPM_PLATFORM.sepchar};

    case 'rootlen'           %-Return length in chars of root directory name
        %=======================================================================
        varargout = {SPM_PLATFORM.rootlen};

    case 'user'                                         %-Return user string
        %=======================================================================
        varargout = {SPM_PLATFORM.user};

    case 'tempdir'                              %-Return temporary directory
        %=======================================================================
        twd = getenv('SPMTMP');
        if isempty(twd)
            twd = tempdir;
        end
        varargout = {twd};


    case {'font','fonts'}    %-Map default font names to platform font names
        %=======================================================================
        if nargin<2, varargout={SPM_PLATFORM.font}; return, end
        switch lower(varargin{2})
            case 'times'
                varargout = {SPM_PLATFORM.font.times};
            case 'courier'
                varargout = {SPM_PLATFORM.font.courier};
            case 'helvetica'
                varargout = {SPM_PLATFORM.font.helvetica};
            case 'symbol'
                varargout = {SPM_PLATFORM.font.symbol};
            otherwise
                warning(['Unknown font ',varargin{2},', using default'])
                varargout = {SPM_PLATFORM.font.helvetica};
        end

    otherwise                                        %-Unknown Action string
        %=======================================================================
        error('Unknown Action string')

        %=======================================================================
end



%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================


function init_platform(comp)             %-Initialise platform variables
%=======================================================================
if nargin<1, comp=computer; end
global SPM_PLATFORM

%-Platform definitions
%-----------------------------------------------------------------------
PDefs = {	'PCWIN',	'win',	0;...
    'MAC2',		'mac',	1;...
    'SUN4',		'unx',	1;...
    'SOL2',		'unx',	1;...
    'HP700',	'unx',	1;...
    'SGI',		'unx',	1;...
    'SGI64',	'unx',	1;...
    'IBM_RS',	'unx',	1;...
    'ALPHA',	'unx',	0;...
    'AXP_VMSG',	'vms',	Inf;...
    'AXP_VMSIEEE',	'vms',	0;...
    'LNX86',	'unx',	0;...
    'GLNX86',	'unx',  0;...
    'VAX_VMSG',	'vms',	Inf;...
    'VAX_VMSD',	'vms',	Inf	};

PDefs = cell2struct(PDefs,{'computer','filesys','endian'},2);


%-Which computer?
%-----------------------------------------------------------------------
ci = find(strcmp({PDefs.computer},comp));
if isempty(ci), error([comp,' not supported architecture for SPM']), end


%-Set bigend
%-----------------------------------------------------------------------
SPM_PLATFORM.bigend = PDefs(ci).endian;
% Commented out as ISIEEE is obsolete and will be removed in future
% versions of MATLAB:
%if ~isieee, SPM_PLATFORM.bigend = Inf; end	%-Last check for IEEE math


%-Set filesys
%-----------------------------------------------------------------------
SPM_PLATFORM.filesys = PDefs(ci).filesys;


%-Set filesystem dependent stuff
%-----------------------------------------------------------------------
%-File separators character
%-Length of root directory strings
%-User name finding
%-(mouse button labels?)
switch (SPM_PLATFORM.filesys)
    case 'unx'
        SPM_PLATFORM.sepchar = '/';
        SPM_PLATFORM.rootlen = 1;
        SPM_PLATFORM.user    = getenv('USER');
    case 'win'
        SPM_PLATFORM.sepchar = '\';
        SPM_PLATFORM.rootlen = 3;
        SPM_PLATFORM.user    = getenv('USERNAME');
        if isempty(SPM_PLATFORM.user)
            SPM_PLATFORM.user = spm_win32utils('username'); end
    case 'mac'
        SPM_PLATFORM.sepchar = ':';
        SPM_PLATFORM.rootlen = 1;			%-** Not sure!?
        SPM_PLATFORM.user    = '';			%-** Dunno!
    otherwise
        error(['Don''t know filesystem ',SPM_PLATFORM.filesys])
end

%-Fonts
%-----------------------------------------------------------------------
switch comp
    case {'SOL2'}	%-Some Sol2 platforms give segmentation violations with Helvetica
        SPM_PLATFORM.font.helvetica = 'Lucida';
        SPM_PLATFORM.font.times     = 'Times';
        SPM_PLATFORM.font.courier   = 'Courier';
        SPM_PLATFORM.font.symbol    = 'Symbol';
    case {'SUN4','SOL2','HP700','SGI','SGI64','IBM_RS','ALPHA','LNX86','GLNX86'}
        SPM_PLATFORM.font.helvetica = 'Helvetica';
        SPM_PLATFORM.font.times     = 'Times';
        SPM_PLATFORM.font.courier   = 'Courier';
        SPM_PLATFORM.font.symbol    = 'Symbol';
    case {'PCWIN'}
        SPM_PLATFORM.font.helvetica = 'Arial Narrow';
        SPM_PLATFORM.font.times     = 'Times New Roman';
        SPM_PLATFORM.font.courier   = 'Courier New';
        SPM_PLATFORM.font.symbol    = 'Symbol';
end



function RT_Version_Callback(hObject, eventdata, handles)
% hObject    handle to RT_Version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RT_Version as text
%        str2double(get(hObject,'String')) returns contents of RT_Version as a double


% --- Executes during object creation, after setting all properties.
function RT_Version_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RT_Version (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function altfolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to altfolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
