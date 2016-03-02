%
% Dartmouth Brain Imaging Center - Conversion Script 
%
% $Id: r2a_convert.m,v 1.2 2006/05/10 15:31:21 jed Exp jed $
%
% Usage: r2convert(parfile,subdir,modality)
%
% valid modality values:
% bold[n], mprage[n], coplanar[n], scout[n]
%

function varargout = r2a_convert(parfile,subdir,modality);
  warning off MATLAB:colon:operandsNotRealScalar;
  
  recfile=parfile;
  if strcmp(recfile(end-2:end),'par')
    recfile(end-2:end)='rec';
  elseif strcmp(recfile(end-2:end),'PAR')
    recfile(end-2:end)='REC';
  end
  
  % create subdirectory if needed
  if exist(subdir) == 0
    if(strcmp(subdir(1),'/') == 1)
      mkdir('/',subdir);
    else
      mkdir(subdir);
    end
  end
  
  Parameters=r2agui('read_par',parfile);
  Precision = strcat ('int', Parameters.bit);
  Precision = char (Precision);
  Vox = Parameters.vox;
  Size = Parameters.dim(1)*Parameters.dim(2)*Parameters.dim(3);
  SizeSlice = Parameters.dim(1)*Parameters.dim(2);
  
  % Open input file
  ID1 = fopen (recfile,'r','l');
  
  Dim = Parameters.dim;
  Type = r2agui('spm_type',Precision);
  Offset = 0;
  Descrip = char (Parameters.name);
  
  Scale = 1;
  Orign = Parameters.fov ./Vox /2 ;
  BytesPerValue=str2num(Parameters.bit)/8;
  
  % Read Data
  iSlice=Parameters.slice_index;
  iSlice(:,12)=iSlice(:,8).*iSlice(:,10).*iSlice(:,11)/8;
  order_slices=iSlice(:,7);
  [os,i]=sort(order_slices); % sort them
  bytespslice_sorted=iSlice(i,12);
  fileposSlice_sorted=[cumsum(bytespslice_sorted)];
  fileposSlice_sorted=[0;fileposSlice_sorted(1:end-1)];
  index_orig=[1:size(order_slices)];
  fileposSlice=fileposSlice_sorted(index_orig(i)); 
  iSlice(:,13)=fileposSlice;
  
  % now sort entire slice_index according to dynamics
  %(fastest varying) and mr_type parameters (slowest varying)
  iSlices_sorted = sortrows(iSlice,[6 3 1]);
  nLine=0;
  dti_increment=0;
  NewFile=[1; (diff(iSlices_sorted(:,3))~=0 | diff(iSlices_sorted( :,6))~=0)]; % determine whether new file has to be opened (change of dynamic/mr_type)
  nr_mrtypes=length(unique(iSlices_sorted(:,6))); % determine number of interleaved image types (e.g. angio)
  while nLine<size(iSlices_sorted,1);
    nLine=nLine+1;
    
    if NewFile(nLine)
      % increment our volume count
      dti_increment=dti_increment+1;
      
      %close previous file
      if nLine>1
        fclose(ID2);
      end
      if nr_mrtypes>1
        mrtype_suffix=num2str(iSlices_sorted(nLine,6),'%04i');
      else
        mrtype_suffix='';
      end
      Precision = strcat ('int', num2str(iSlices_sorted(nLine,8)));
      Precision = char (Precision);
      cDim=[iSlices_sorted(nLine,10:11),Dim(3)];
      % Test for single volume
      
      if Parameters.dyn == 1
        VolName =  [subdir '/' modality mrtype_suffix '.img'];
      else
        
        if modality(1:3) == 'dti'
          VolName =  [subdir '/' modality mrtype_suffix '_' ...
                      num2str(dti_increment,'%04i') '.img'];
        else
          VolName =  [subdir '/' modality mrtype_suffix '_' ...
                      num2str(iSlices_sorted(nLine,3),'%04i') '.img'];
        end
      end
      
      ID2 = fopen (VolName, 'w');
      P = VolName;
      r2agui('spm_hwrite',P,cDim,Vox,Scale,Type,Offset,round(Orign),Descrip);
      % Test for single volume
      if Parameters.dyn == 1
        disp (['Writing file: ' subdir '/' modality]);
      else
        disp (['Writing file: ' subdir '/' modality,' volume:',num2str(iSlices_sorted(nLine,3),'%04.0f')]);
      end
    end
    
    fseek(ID1,iSlices_sorted(nLine,13),-1);
    SliceData = fread (ID1,cDim(1)*cDim(2), Precision);
    fwrite(ID2, SliceData,Precision);
  end
  fclose(ID2);
  fclose(ID1);
  
