%
% Dartmouth Brain Imaging Center
%
% $Id: par_boil.m,v 1.3 2008/09/12 18:39:41 jed Exp jed $
%
% Usage: par_boil (indir)
% 

function varargout = par_boil(indir);

% get working directory
cwd=pwd;

if nargin == 0
 disp('Need input and output directories!');
 return
end

% check for existence of the input directory
if exist(indir) ~= 0
  cd(indir)
 else
  disp('Input directory does not exist!');
  return;
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
 
  % obtain modality and parameters
  [modality,Parameters] = par_func(parfile); 

  name = char(Parameters.name);

  switch lower(modality)
   case 'scout'
    if scout_i == 1
     modality='scout';
    else
     modality=strcat('scout',num2str(scout_i));
    end
    scout_i=scout_i + 1;
   case 'dti'
    if dti_i == 1
     modality='dti';
    else
     modality='dti'; %strcat('dti',num2str(dti_i));
    end
    dti_i=dti_i + 1;
   case 'mprage'
    if mprage_i == 1
     modality='mprage';
    else
     modality=strcat('mprage',num2str(mprage_i));
    end
   mprage_i = mprage_i + 1;
   case 'bold'
    if bold_i == 1
     modality='bold1';
    else
     modality=strcat('bold',num2str(bold_i));
    end
   bold_i=bold_i + 1;
  case 'coplanar'
    if coplanar_i == 1
     modality='coplanar';
    else
     modality=strcat('coplanar',num2str(coplanar_i));
    end
   coplanar_i=coplanar_i + 1;
  end

  output=sprintf('Scan: %s  Modality: %s (%s)',parfile,modality,name);
  disp(output);

  i=i+1;
end

% change back
cd(cwd);
