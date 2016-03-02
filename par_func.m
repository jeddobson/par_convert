%
% Dartmouth Brain Imaging Center
%
% $Id: par_func.m,v 1.2 2009/01/14 20:14:45 jed Exp jed $
%
% Usage: par_func(object)
% 
% returns modality of PAR/REC file
%
%  Modified by Yu-CHien Wu 2009/07/21
%    change the logic of assigning modality name 
%    based on technique name, dynamics and diffusion flag 
%    (instead of slice and volume number)  
%

function [modality,Parameters] = par_func(parfile);


% get working directory
cwd=pwd;

Parameters=r2agui('read_par',parfile);
volumes = Parameters.dyn;
name = char(Parameters.name);
slices = Parameters.slice;
voxel_sz = Parameters.vox(1); % look just at dim1, silly?

% moved up for older data without additional fields (Jed 02/07/12)
if strcmp(name,'Sense Ref')
   modality='sense';
   return;
end

% Wu added
tech = Parameters.tech;

Diffusion = Parameters.diffusion;
bvalno = Parameters.bvalno;
bvecno = Parameters.bvecno;

% For all not specified in the following 
% name after MR technique: researchers should know what they are doing
modality=tech;

% check for Scout
%if strcmp(name,'Scout')
% 03/11/2011: Added Survey for 32ch protocol (Jed)
if isempty(strfind(name,'Scout'))~=1 || isempty(strfind(name,'scout'))~=1 || isempty(strfind(name,'Survey'))~=1
   modality='scout';
   return;
end

if strcmp(name,'Sense Ref')
   modality='sense';
   return;
end

% COPLANAR
% coplanar is not a standard MR sequence and could use any kinds of technique, voxel resolution, slice number 
% it usually accompany with other functional study
% there is only one group use coplanar as far as I know (Catherine Norris)
% 
if isempty(strfind(name,'Coplanar'))~=1 || isempty(strfind(name,'coplanar'))~=1
  modality='coplanar';
  return;
end

% DTI, many volumes, smaller voxels -- DTI before MPRAGE
if Diffusion == '1' 
    if isempty(bvalno)~=1   % added by Jed 02/15/2011 to convert old dti data
      if bvalno==2 && bvecno >=6
        modality='dti';
      else
        modality='dwi';
      end
    else
     % pre-V4.X
    modality='dti';
    end
  return;
end

% 32channel hires

% MPRAGE: High resolution images: meanc voxel size is small    ..Wu 
% T1W, single volume, higher slices
if strcmp(tech,'T1TFE') && voxel_sz < 1.5 
  modality='mprage';
  return;
end
  
% fMRI, many volumes, bigger voxels
if strcmp(tech,'FEEPI') && volumes > 1 
    modality='bold';
  return;
end


return;

