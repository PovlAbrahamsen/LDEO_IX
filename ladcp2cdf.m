function [] = ladcp2cdf(fname,dr_struct,da,p,ps,f,att);
% function [] = ladcp2cdf(fname,dr_struct,da,p,ps,f,att);
%
% function to save LADCP data into a netcdf file for MatLab version 2012a
%
% input  :	fname		- output filename
%			dr_struct	- main inversion results (velocity profiles)
%			da..att	- arbitrary metadata structures
%
% Created By:   Diana Cardoso, Bedford Institute of Oceangraphy
%               Diana.Cardoso@dfo-mpo.gc.ca
% Description:  Based on LDEO software to Process LADCP, version IX.8,
%               script ladcp2cdf.m version 0.1	last change 08.03.2002. 
%               maintained by A.M. Thurnherr and downloaded from:
%       http://www.ldeo.columbia.edu/cgi-bin/ladcp-cgi-bin/hgwebdir.cgi
%       The function ladcp2cdf was changed to run with the the Matlab
%       version 2012, which now supports netcdf.

%======================================================================
%                    L A D C P 2 C D F . M 
%                    doc: Thu Aug 15 10:52:55 2013
%                    dlm: Wed Aug 28 12:31:16 2013
%                    (c) 2013 A.M. Thurnherr, from code contributed by D. Cardoso
%                    uE-Info: 99 0 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% NOTES:
%	- This version creates slightly different files than the original version
%	  created by Visbeck/Krahmann. In the original version, the contents of the
%	  dr structure end up as top-level variables and the contents of
%	  the da, p, ps, f and att structures end of as global attributes. In 
%	  the new version, the latter are saved as sub-structures, with _struct appended
%	  to the internal names to avoid conflicts.

% HISTORY:
%   Aug 15, 2013: - incorporated this code, supplied Diana Cardoso, into IX_10beta
%		  - modified doc in header
%		  - renamded struct variable to dr_struct
%		  - removed 'cd' in and out of results directory (pathnames work just fine)
%		  - delete netcdfile before it is written to (old 'clobber' option)
%		  - removed 'l' suffix from all dims
%		  - replaced yes/no logical vals by true/false
%		  - renamed substructures from st2..st6 to internal names (da,p,ps,f,att)
%   Aug 28, 2013: - incorporated bug fix provided by Diana Cardoso to prevent lat,lon,name and 
%		    date to be stored 2cd in the nc file, which can make the code
%		    bomb if the length of any other var is 6 (or equal to the length of name?)

% check arguments
if nargin<2
  error('need two input arguments')
end
if ~isstruct(dr_struct)
  error('second argument must be a dr structure')
end

netcdfile = deblank(fname); %remove any blanks from string end
if exist(netcdfile,'file')
	delete(netcdfile)
end

%Determine dimensions of variables
lbot = 0;
lz = 0;
ltim = 0;
lsadcp = 0;

if isfield(dr_struct,'z');
  lz = length(getfield(dr_struct,'z'));
end  
if isfield(dr_struct,'tim');
  ltim = length(getfield(dr_struct,'tim'));
end  
if isfield(dr_struct,'zbot');
  lbot = length(getfield(dr_struct,'zbot'));
end  
if isfield(dr_struct,'z_sadcp');
  lsadcp = length(getfield(dr_struct,'z_sadcp'));
end  

% % define dimensions in netcdf file and standard variables
nccreate(netcdfile,'lat','Dimensions',{'lat' 1},'Datatype','single');
nccreate(netcdfile,'lon','Dimensions',{'lon' 1},'Datatype','single');
nccreate(netcdfile,'date','Dimensions',{'date' 6},'Datatype','int32');
nccreate(netcdfile,'name','Dimensions',{'name' length(getfield(dr_struct,'name'))},'Datatype','char');

% store standard variables
ncwrite(netcdfile,'lat',dr_struct.lat);
ncwrite(netcdfile,'lon',dr_struct.lon);
ncwrite(netcdfile,'date',dr_struct.date);
ncwrite(netcdfile,'name',dr_struct.name);

% parse fieldnames, define the proper variable and store it
fnames = fieldnames(dr_struct);
nn=strncmp('name',fnames,6); nla=strncmp('lat',fnames,3); 	% find, name, lat,lon,date from fnames
nlo=strncmp('lon',fnames,3); nda=strncmp('date',fnames,4);
ntot=[nn+nla+nlo+nda]; Ktot = logical(ntot);
fnames(Ktot,:)=[];						% remove , name, lat,lon,date from fnames

for n=1:size(fnames,1)
  dummy = getfield(dr_struct,fnames{n});
  if length(dummy)==lz
    nccreate(netcdfile,fnames{n},'Dimensions',{fnames{n} fix(lz)},'Datatype','single');
    ncwrite(netcdfile,fnames{n},dummy);
  end
  if length(dummy)==ltim
    nccreate(netcdfile,fnames{n},'Dimensions',{fnames{n} fix(ltim)},'Datatype','single');
    ncwrite(netcdfile,fnames{n},dummy);
  end
  if length(dummy)==lbot
    nccreate(netcdfile,fnames{n},'Dimensions',{fnames{n} fix(lbot)},'Datatype','single');
    ncwrite(netcdfile,fnames{n},dummy);
  end
  if length(dummy)==lsadcp
    nccreate(netcdfile,fnames{n},'Dimensions',{fnames{n} fix(lsadcp)},'Datatype','single');
    ncwrite(netcdfile,fnames{n},dummy);
  end
end

add_struct(netcdfile,'da_struct',da)
add_struct(netcdfile,'p_struct',p)
add_struct(netcdfile,'ps_struct',ps)
add_struct(netcdfile,'f_struct',f)
add_struct(netcdfile,'att_struct',att)

end % function

%----------------------------------------------------------------------

function [] = add_struct(ncf,snm,a)
  
   fnames = fieldnames(a);
   if isstruct(a)
      if ~isstruct(eval(['a.' fnames{1}])) % No SubStructure
	nccreate(ncf,snm,'Datatype','char');

	for n = 1:size(fnames,1)
		dummy = getfield(a,fnames{n});
	    	if size(dummy,1)==1
	       		if isstr(dummy)
                		ncwriteatt(ncf,snm,fnames{n},dummy);
	       		elseif islogical(dummy)
		        	if dummy, dummy='true';
		               	else, 	  dummy='false';
		               	end
				ncwriteatt(ncf,snm,fnames{n},dummy);
			else
		                ncwriteatt(ncf,snm,fnames{n},dummy(:));
			end    
		end
	end % for n
      else % SubStructures -> Variable Attributes
	for n = 1:size(fnames,1)
		atts = eval(['fieldnames(a.' fnames{n} ');']);
	        finfo = ncinfo(ncf);
        	FieldNames = {finfo. Variables.Name};
	        existField=strmatch(fnames{n}, FieldNames);
        	if isempty(existField)
			nccreate(ncf,fnames{n},'Datatype','char');
	        end
		for j = 1:size(atts,1)
			dummy = eval(['a.' fnames{n} '.' atts{j} ';']);
			if size(dummy,1) == 1
				if ischar(dummy)
		                    ncwriteatt(ncf,fnames{n},atts{j},dummy);
                		else
		                    ncwriteatt(ncf,fnames{n},atts{j},dummy(:));
                		end
			end
		end	       
	 end % for n
      end % else (substructures or not)
   else % if issstruct(a)
      disp(' not structure')
   end
end % function

