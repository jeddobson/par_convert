%
% Dartmouth Brain Imaging Center
%
% $Id: par_convert_afni.m,v 1.2 2008/03/13 19:13:32 jed Exp jed $
%
% Usage: par_convert(indir,outdir)
% 
% This script will convert a directory of PAR/REC files to AFNI.
% You can manually convert a single PAR/REC file pair with r2a_convert
% or by using r2agui directly (the script called by r2a_convert).
%


function varargout = par_convert(indir,outdir);

disp('starting afni conversion...');

% get working directory
cwd=pwd;

if nargin == 0
 disp('Need input and output directories!');
 return
end

% check for existence of the input directory
if exist(indir) ~= 0
  cd(indir);
 else
  error=sprintf('Input directory (%s) does not exist!',indir);
  disp(error);
  return;
end

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

% Philips DTI files have lower case?
parfiles = cat(1,parfiles1,parfiles2);


i=1;
while i<size(parfiles,1)+1;
  parfile=parfiles(i).name;
  
  %
  % Replicate my Perl logic for figuring what is what...
  %
  
  Parameters=r2agui('read_par',parfile);
  volumes = Parameters.dyn;
  slices = Parameters.slice;
  voxel_sz = Parameters.vox(1); % look just at dim1, silly?
  
  % Scout, single volume, fewer slices
  if volumes == 1 && slices < 15
    if scout_i == 1
     modality='scout';
    else
     modality=strcat('scout',num2str(scout_i));
    end
    subdir='ANATOMY';
    subdir=strcat(outdir,'/',subdir);
    r2a_convert(parfile,subdir,modality);
    scout_i=scout_i+1;
  end
  
  % MPRAGE, single volume, higher slices
  if volumes == 1 && slices > 90
    if mprage_i == 1
    modality='mprage';
    else
    modality=strcat('mprage',num2str(mprage_i));
    end
    subdir='ANATOMY';
    subdir=strcat(outdir,'/',subdir);
    r2a_convert(parfile,subdir,modality);

    % reorient - xz, yz
    cmd=sprintf('reorient %s/%s %s/%s xz o > /dev/null 2>&1',subdir, ...
                modality,subdir, modality);
    unix(cmd);
    
    cmd=sprintf('reorient %s/%s %s/%s yz o > /dev/null 2>&1',subdir, ...
                modality,subdir, modality);
    unix(cmd);
    
    wd = pwd;
    cd(subdir);
    cmd=sprintf('to3d -anat -orient LPI -session . -prefix %s %s.hdr', ...
                modality,modality);
    unix(cmd);
    cmd=sprintf('rm %s*.{hdr,img}',modality);
    unix(cmd);
    cd(wd);

    mprage_i=mprage_i+1;
  end
  
  % COPLANAR, single volume, few slices.
  if volumes == 1 && slices > 15 & slices < 90
    if coplanar_i == 1 
       modality='coplanar';
   else
       modality=strcat('coplanar',num2str(coplanar_i));
    end
    subdir='ANATOMY';
    subdir=strcat(outdir,'/',subdir);
    r2a_convert(parfile,subdir,modality);
    coplanar_i=coplanar_i+1;
  end
  
  % fMRI, many volumes, bigger voxels
  if volumes > 1 && voxel_sz > 2
    modality=strcat('bold',num2str(bold_i));
    subdir='FUNCTIONAL';
    subdir=strcat(outdir,'/',subdir);
    r2a_convert(parfile,subdir,modality);
    wd = pwd;
    cd(subdir);
    % reorient y

    cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file y o > /dev/null 2>&1;done',...
    subdir,modality);
    unix(cmd);

    % reorient x
    cmd=sprintf('for file in %s/%s*.hdr;do reorient $file $file x o > /dev/null 2>&1;done',...
    subdir,modality);
    unix(cmd);

    cmd=sprintf('to3d -time:zt %d %d 2012 alt+z -orient LPI -session . -prefix %s  %s*.hdr', ...
                slices,volumes,modality,modality);
    unix(cmd);
    cmd=sprintf('rm %s*.{hdr,img}',modality);
    unix(cmd);
    cd(wd);
    bold_i=bold_i+1;
  end

  % DTI, many volumes, smaller voxels
  if volumes > 1 && voxel_sz < 2
    modality=strcat('dti',num2str(dti_i));
    subdir='DTI';
    subdir=strcat(outdir,'/',subdir);
    r2a_convert(parfile,subdir,modality);
    dti_i=dti_i+1;
  end
  
  i=i+1;
end

% change back
cd(cwd);
