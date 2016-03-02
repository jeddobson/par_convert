%
% Dartmouth Brain Imaging Center
%
% $Id: par_convert.m,v 1.26 2014/12/22 19:41:51 jed Exp jed $
% Modified by Yu-Chien Wu 2009/07/20
%
% This script will convert a directory of PAR/REC files to SPM/ANALYZE.
% You can manually convert a single PAR/REC file pair with r2a_convert
% or by using r2agui directly (the script called by r2a_convert).
%
% Supported Formats: nifti, analyze, afni
%

function par_convert(indir,outdir,format)
version = '$Revision: 1.26 $';
[rc,sysuser]=system('whoami');

if nargin == 0
    disp('Need input and output directories!');
    return
end

% If format is not specified, assume analyze
if nargin < 3
    format = 'analyze';
end

% get working directory
cwd=pwd;

% check to see if output is relative or absolute
% and convert to absolute
if outdir(1) ~= '/'
 outdir=strcat(cwd,'/',outdir);
end

% check for existence of the input directory
if exist(indir) ~= 0
    cd(indir)
else
    disp('Input directory does not exist!');
    return;
end

% check for existence of the output directory
outdir_preproc=strcat(outdir,'/SCRIPTS');
if exist(outdir) == 7 && exist(outdir_preproc) == 7
    disp('Output directory already exists: using pre-proc--okay');
elseif exist(outdir) == 7 && exist(outdir_preproc) == 0
    disp('Output directory already exists');
    cd(cwd); 
    return;
end
    

% Create cell for bold sizes
boldarray=cell(10,2);

% Display selected output format
o=sprintf('Writing data as: %s',format);
disp(o);

% Setup logging
message=sprintf('logger -i "par_convert: %s/%s"',indir,format);
system(message);
logdir=strcat(outdir,'/LOG');
logfile=strcat(logdir,'/','par_convert.log');
td=date;
mkdir('/',logdir);
fid = fopen(logfile,'w');
fprintf(fid,'Using par_convert version: %s\n',version);
fprintf(fid,'Date: %s\n',td);
fprintf(fid,'User: %s\n\n',sysuser);
fprintf(fid,'Output format: %s\n',format);
fprintf(fid,'Input Directory: %s/%s\n',cwd,indir);
fprintf(fid,'Output Directory: %s\n',outdir);

% create subdirectory if needed
if exist(outdir) == 0
    mkdir('/',outdir);
end

% get a list of all PAR files
parfiles1 = dir('*.PAR');
parfiles2 = dir('*.par');

% set indices to one
scout_i=1;
mprage_i=1;
coplanar_i=1;
bold_i=1;
dti_i=1;
dwi_i=1;

% Philips DTI files have lower case?
parfiles = cat(1,parfiles1,parfiles2);

i=1;
while i<size(parfiles,1)+1;
    parfile=parfiles(i).name;
    
    
    %
    % Replicate my Perl logic for figuring what is what...
    %
    % Determine the run number and sort into proper directories defined by
    %
    
    % obtain modality and parameters
    [modality,Parameters] = par_func(parfile);     %% Wu modify
    par_version = char(Parameters.ResToolsVersion);
    name = char(Parameters.name);
    volumes = char(Parameters.dyn);
    slices = char(Parameters.slice);
    
    
    % write to log
    fprintf(fid,'\n---%s---\n',parfile);
    
    % log the size
    [x,y] = size(parfile);
    recfile=strcat(parfile(1: y-4),'.REC');
    fd=dir(recfile);
    
    
    fprintf(fid,'date: %s\n',fd.date);
    fprintf(fid,'size: %d\n\n',fd.bytes);
    fprintf(fid,'Protocol Name: %s\n',name);
    fprintf(fid,'Volumes: %d\n',volumes);
    fprintf(fid,'Slices: %d\n',slices);
    
    switch lower(modality)
        
        case 'sense'
            modality='sense'; % don't convert
            
            % scout
        case 'scout'
            if scout_i == 1
                modality='scout';
            else
                modality=strcat('scout',num2str(scout_i));
            end
            subdir='ANATOMY';
            subdir=strcat(outdir,'/',subdir);
            fprintf(fid,'Converting as: %s\n',modality);
            r2a_convert(parfile,subdir,modality);
            scout_i=scout_i + 1;
            
            % mprage
        case 'mprage'
            if mprage_i == 1
                modality='mprage';
            else
                modality=strcat('mprage',num2str(mprage_i));
            end
            subdir='ANATOMY';
            subdir=strcat(outdir,'/',subdir);

           %log 
           fprintf(fid,'Converting as: %s (%s)\n',modality,Parameters.sliceorient);
           r2a_convert(parfile,subdir,modality);

	   if isempty(strfind(Parameters.sliceorient,'2'))~=1
            % reorient - xz, yz
            cmd=sprintf('reorient %s/%s %s/%s xz o > /dev/null 2>&1',subdir, ...
                modality,subdir, modality);
            unix(cmd);
            
            cmd=sprintf('reorient %s/%s %s/%s yz o > /dev/null 2>&1',subdir, ...
                modality,subdir, modality);
            unix(cmd);

           else % rotation around z axis (was mirror) for sliceorient == 1
            %reorient - zz
	    fprintf(fid,'rotating mprage');
            cmd=sprintf('reorient %s/%s %s/%s zz o > /dev/null 2>&1',subdir, ...
              modality,subdir, modality);
            unix(cmd);

          end 

            mprage_i = mprage_i + 1;
            
            % format changes
            if strcmp(format,'afni')
                wd=pwd;
                cd(subdir);
                cmd=sprintf('to3d -anat -orient LPI -session . -prefix %s %s.hdr 2>&1 | tee -a %s/afni_convert.log', ...
                    modality,modality,logdir);
                unix(cmd);
                cmd=sprintf('rm %s*.{hdr,img}',modality);
                unix(cmd);
                cd(wd);
            end
            
            if strcmp(format,'nifti') || strcmp(format,'nifti3d') 
                wd=pwd;
                cd(subdir);
                
                % convert to NIfTI-1
                cmd=sprintf('fslchfiletype NIFTI_GZ %s.nii.gz %s.hdr', ...
                    modality,modality);
                unix(cmd);
                cmd=sprintf('rm %s*.{hdr,img}',modality);
                unix(cmd);
                
                % force qcode & scode to be set to 1
                % otherwise qform & sform are not used)
                cmd=sprintf('fslorient -setsformcode 1 %s.nii.gz; fslorient -setqformcode 1 %s.nii.gz', modality,modality);
                unix(cmd);
                
                % force header to record neurological orientation
                cmd=sprintf('fslorient -forceneurological %s.nii.gz', modality);
                unix(cmd);
                
                cd(wd);
            end
            
        case 'coplanar'
            if coplanar_i == 1
                modality='coplanar';
            else
                modality=strcat('coplanar',num2str(coplanar_i));
            end
            subdir='ANATOMY';
            subdir=strcat(outdir,'/',subdir);
            % log
            fprintf(fid,'Converting as: %s\n',modality);
            r2a_convert(parfile,subdir,modality);
            
            % reorient - mirror y
            cmd=sprintf('reorient %s/%s %s/%s y o > /dev/null 2>&1',subdir,modality,subdir,modality);
            unix(cmd);
            cmd=sprintf('reorient %s/%s %s/%s x o > /dev/null 2>&1',subdir,modality,subdir,modality);
            unix(cmd);
            
            % format
            if strcmp(format,'nifti') || strcmp(format,'nifti3d')
                wd=pwd;
                cd(subdir);
                
                cmd=sprintf('fslchfiletype NIFTI_GZ %s.nii.gz %s.hdr', ...
                    modality,modality);
                unix(cmd);
                cmd=sprintf('rm %s*.{hdr,img}',modality);
                unix(cmd);
                
                % force qcode & scode to be set to 1
                % otherwise qform & sform are not used)
                cmd=sprintf('fslorient -setsformcode 1 %s.nii.gz; fslorient -setqformcode 1 %s.nii.gz', modality,modality);
                unix(cmd);
                
                % force header to record neurological orientation
                cmd=sprintf('fslorient -forceneurological %s.nii.gz', modality);
                unix(cmd);
                
                
                cd(wd);
            end
            if strcmp(format,'afni')
                wd=pwd;
                cd(subdir);
                cmd=sprintf('to3d -anat -orient LPI -session . -prefix %s %s.hdr 2>&1 | tee -a %s/afni_convert.log', ...
                    modality,modality,logdir);
                unix(cmd);
                cmd=sprintf('rm %s*.{hdr,img}',modality);
                unix(cmd);
                cd(wd);
            end
            coplanar_i=coplanar_i + 1;
            
            % fMRI
        case 'bold'
            if bold_i == 1
                modality='bold1';
            else
                modality=strcat('bold',num2str(bold_i));
            end
            boldarray{bold_i,2} = fd.bytes;
            boldarray{bold_i,1} = modality;
            subdir='FUNCTIONAL';
            subdir=strcat(outdir,'/',subdir);
            
            % log
            fprintf(fid,'Converting as: %s\n',modality);
            r2a_convert(parfile,subdir,modality);
            
            % reorient y
            cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file y o ;done',subdir,modality);
            unix(cmd);
            
            % reorient x
            cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file x o ;done',subdir,modality);
            unix(cmd);
            
            % format
            if strcmp(format,'afni')
                wd=pwd;
                cd(subdir);
                tr=Parameters.RT*1000;
                cmd=sprintf('/afs/dbic.dartmouth.edu/usr/pkg/afni/AFNI_2011_12_21/amd64_linux26/bin/to3d -time:zt %d %d %d alt+z -orient LPI -session . -prefix %s %s*.hdr 2>&1 | tee -a %s/afni_convert.log', ...
                    slices,volumes,tr,modality,modality,logdir);
		fprintf(fid,'AFNI command: %s\n',cmd);	
                unix(cmd);
                cmd=sprintf('rm %s*.{hdr,img}',modality);
                unix(cmd);
                cd(wd);
            end
            
            if strcmp(format,'nifti') || strcmp(format,'nifti3d')
                wd=pwd;
                cd(subdir);
                
                cmd=sprintf('FSLOUTPUTTYPE=NIFTI fslmerge -t %s.nii  %s*.hdr', ...
                    modality,modality);
                unix(cmd);
                
                % force qcode & scode to be set to 1
                % otherwise qform & sform are not used)
                cmd=sprintf('fslorient -setsformcode 1 %s.nii; fslorient -setqformcode 1 %s.nii', modality,modality);
                unix(cmd);
                
                % force header to record neurological orientation
                cmd=sprintf('fslorient -forceneurological %s.nii', modality);
                unix(cmd);
               
		if strcmp(format,'nifti') 
                  % add tr + slice timing correction & compress
                  cmd=sprintf('fixniftihdr %s/%s %d',subdir,modality,Parameters.RT);
                  unix(cmd);
                else
                  % add tr + slice timing correction & compress
                  cmd=sprintf('fixniftihdr3d %s/%s %d',subdir,modality,Parameters.RT);
                  unix(cmd);
		end
                
                cmd=sprintf('rm %s*.{hdr,img}',modality);
                unix(cmd);
                cd(wd);
            end


            bold_i=bold_i + 1;
            
            % DTI
        case 'dti'
            if dti_i == 1
                modality='dti1'
            else
                %modality='dti';
                modality=strcat('dti',num2str(dti_i))
            end
            
            subdir=['DTI' num2str(dti_i)];
            subdir=strcat(outdir,'/',subdir);
            
            % log
            fprintf(fid,'Converting as: %s\n',modality);
            
            if strcmp(par_version,'V4')
                
                % convert
                r2a_convert(parfile,subdir,modality);
                
                % reorient y
                cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file y o > /dev/null 2 >&1;done',...
                    subdir,modality);
                unix(cmd);
                
                % reorient x
                cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file x o > /dev/null 2 >&1;done',...
                    subdir,modality);
                unix(cmd);
                
            else
                filelist={parfile};
                options.prefix=['dti' num2str(dti_i)];
                options.usealtfolder=1;
                options.pathpar=strcat(pwd,'/');
                options.subaan=1;
                options.altfolder=outdir;
                   
                options.angulation=1;
                options.rescale=1;
                options.usefullprefix=1;
                options.outputformat=2;
                
                output=convert_r2a(filelist,options);
                
                % reorient y
                cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file y o >/dev/null 2>&1 ;done',...
                    subdir,modality);
                unix(cmd);
                
                % reorient x
                cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file x o >/dev/null 2>&1;done',...
                    subdir,modality);
                unix(cmd);
                
                
                %       % create output directory by hand (sort of)
                %      mkdir('/',subdir);
                %
                %       cmd=sprintf('3dPAR2AFNI -n -o %s %s > /dev/null 2>&1',subdir,parfile);
                %       unix(cmd);
                %       root=parfile(1:end-4);
                %
                %       if strcmp(format,'analyze')
                %         cmd=sprintf('FSLOUTPUTTYPE=NIFTI_PAIR fslsplit %s/%s.nii %s/%s_',subdir, ...
                %                     root,subdir,modality);
                %         unix(cmd);
                %
                %	% reorient data
                %        cmd=sprintf('for file in %s/%s*.hdr;do fslswapdim $file -x -y z $file; done',...
                %		subdir,modality);
                %        unix(cmd);
                
                %	cmd=sprintf('rm %s/%s*.{img,hdr}',subdir,modality);
                %	unix(cmd);
                
                % convert to analyze
                %        cmd=sprintf('for file in %s/%s*.gz;do fslchfiletype ANALYZE $file > /dev/null; done',...
                %                     subdir,modality);
                %        unix(cmd);
                %
                %
                %       elseif strcmp(format,'nifti') % with nifti the orientation is going to
                %                                     % be off a bit until we sort out how to
                %                                     % store qform & sform data correctly.
                %         cmd=sprintf('FSLOUTPUTTYPE=NIFTI_GZ fslsplit %s/%s.nii %s/%s_',subdir, ...
                %                     root,subdir,modality);
                %         unix(cmd);
                %       end
                %
                %       cmd=sprintf('rm %s/%s*',subdir,root);
                %       unix(cmd);
            end
            dti_i = dti_i + 1;
            
            % Modified by Yu-Chien Wu July 27, 2009
            % DWI: all other diffusion weighted scans, except diffusion
            % tensor
        case 'dwi'
            if dwi_i == 1
                modality='dwi1';
            else
                modality=strcat('dwi',num2str(dwi_i));
            end
            
            subdir=['DWI' num2str(dwi_i)];
            subdir=strcat(outdir,'/',subdir);
            
            % log
            fprintf(fid,'Converting as: %s\n',modality);
            
            if strcmp(par_version,'V4')
                
                % convert
                r2a_convert(parfile,subdir,modality);
                
                % reorient y
                cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file y o > /dev/null 2 >&1;done',...
                    subdir,modality);
                unix(cmd);
                
                % reorient x
                cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file x o > /dev/null 2 >&1;done',...
                    subdir,modality);
                unix(cmd);
                
            else
                filelist={parfile};
                options.prefix=['dwi' num2str(dwi_i)];
                
                options.altfolder=outdir;
                options.subaan=1;
                options.usealtfolder=1;
                              
                options.pathpar=strcat(pwd,'/');
                options.angulation=1;
                options.rescale=1;
                options.usefullprefix=1;
                options.outputformat=2;
                output=convert_r2a(filelist,options);
                
                % reorient y
                cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file y o >/dev/null 2>&1 ;done',...
                    subdir,modality);
                unix(cmd);
                
                % reorient x
                cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file x o >/dev/null 2>&1;done',...
                    subdir,modality);
                unix(cmd);
                
                
                %       % create output directory by hand (sort of)
                %      mkdir('/',subdir);
                %
                %       cmd=sprintf('3dPAR2AFNI -n -o %s %s > /dev/null 2>&1',subdir,parfile);
                %       unix(cmd);
                %       root=parfile(1:end-4);
                %
                %       if strcmp(format,'analyze')
                %         cmd=sprintf('FSLOUTPUTTYPE=NIFTI_PAIR fslsplit %s/%s.nii %s/%s_',subdir, ...
                %                     root,subdir,modality);
                %         unix(cmd);
                %
                %	% reorient data
                %        cmd=sprintf('for file in %s/%s*.hdr;do fslswapdim $file -x -y z $file; done',...
                %		subdir,modality);
                %        unix(cmd);
                
                %	cmd=sprintf('rm %s/%s*.{img,hdr}',subdir,modality);
                %	unix(cmd);
                
                % convert to analyze
                %        cmd=sprintf('for file in %s/%s*.gz;do fslchfiletype ANALYZE $file > /dev/null; done',...
                %                     subdir,modality);
                %        unix(cmd);
                %
                %
                %       elseif strcmp(format,'nifti') % with nifti the orientation is going to
                %                                     % be off a bit until we sort out how to
                %                                     % store qform & sform data correctly.
                %         cmd=sprintf('FSLOUTPUTTYPE=NIFTI_GZ fslsplit %s/%s.nii %s/%s_',subdir, ...
                %                     root,subdir,modality);
                %         unix(cmd);
                %       end
                %
                %       cmd=sprintf('rm %s/%s*',subdir,root);
                %       unix(cmd);
            end
            dwi_i = dwi_i + 1;
            
            
            % other MR methods: wu added July 27, 2009
        otherwise
            % classify by series number  
            subdir=['S_' num2str(i)];
            subdir=strcat(outdir,'/',subdir);
            % log
            fprintf(fid,'Converting as: %s\n',modality);
            r2a_convert(parfile,subdir,modality);
            
            % reorient - xz, yz
            cmd=sprintf('reorient %s/%s %s/%s xz o > /dev/null 2>&1',subdir, ...
                modality,subdir, modality);
            unix(cmd);
            
            cmd=sprintf('reorient %s/%s %s/%s yz o > /dev/null 2>&1',subdir, ...
                modality,subdir, modality);
            unix(cmd);
            
            
%             % format changes
%             if strcmp(format,'afni')
%                 wd=pwd;
%                 cd(subdir);
%                 cmd=sprintf('to3d -anat -orient LPI -session . -prefix %s %s.hdr 2>&1 | tee %s/afni_convert.log', ...
%                     modality,modality,logdir);
%                 unix(cmd);
%                 cmd=sprintf('rm %s*.{hdr,img}',modality);
%                 unix(cmd);
%                 cd(wd);
%             end
%             
%             if strcmp(format,'nifti')
%                 wd=pwd;
%                 cd(subdir);
%                 
%                 % convert to NIfTI-1
%                 cmd=sprintf('fslchfiletype NIFTI_GZ %s.nii.gz %s.hdr', ...
%                     modality,modality);
%                 unix(cmd);
%                 cmd=sprintf('rm %s*.{hdr,img}',modality);
%                 unix(cmd);
%                 
%                 % force qcode & scode to be set to 1
%                 % otherwise qform & sform are not used)
%                 cmd=sprintf('fslorient -setsformcode 1 %s.nii.gz; fslorient -setqformcode 1 %s.nii.gz', modality,modality);
%                 unix(cmd);
%                 
%                 % force header to record neurological orientation
%                 cmd=sprintf('fslorient -forceneurological %s.nii.gz', modality);
%                 unix(cmd);
%                 
%                 cd(wd);
%             end
            
    % -------- end case ----------------%        
            
    end
    
    % final increment & end
    i=i+1;
end

% do a final bold size check
b=cat(2,boldarray{:,2});
if std(b) ~= 0 && isnan(std(b)) ~= 1
    disp('Potential problem with fMRI runs, please check log.');
    boldarray
    fprintf(fid,'\n\nERROR: fMRI sizes do not match\n');
end


% close logfile
fclose(fid);

%change back
cd(cwd);


