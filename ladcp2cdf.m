%======================================================================
%                    L A D C P 2 C D F _ V 2 . M 
%                    doc: Fri Jul 24 09:31:43 2015
%                    dlm: Mon Jul 27 15:47:29 2015
%                    (c) 2015 A.M. Thurnherr
%                    uE-Info: 21 0 NIL 0 0 72 2 2 8 NIL ofnI
%======================================================================

% HISTORY:
%  Jul 24, 2015: - created, losely based on code provided by D. Cardoso
%  Jul 27, 2015: - reverted to Martin's original dimensions as requested
%		   by E. Firing

function [] = ladcp2cdf(fname,dr,da,p,ps,f,att);
    netcdfile = deblank(fname);
    if exist(netcdfile,'file'), delete(netcdfile); end

    add_dr_struct(netcdfile,'dr',dr,att)
    add_struct(netcdfile,'da',da)
    add_struct(netcdfile,'p',p)
    add_struct(netcdfile,'ps',ps)
    add_struct(netcdfile,'f',f)
end

%----------------------------------------------------------------------

function [] = addatts(ncf,fnm,att)
    as = getfield(att,fnm);
    an = fieldnames(as);
    for i=1:length(an)
    	ncwriteatt(ncf,sprintf('dr.%s',fnm),char(an(i)),char(getfield(as,char(an(i)))));
    end
end   

function [] = add_dr_struct(ncf,snm,dr,att)

    % scalars & misc dims
    nccreate(ncf,'dr.name','Datatype','char','Dimensions',{'dr.name' length(dr.name)});
	ncwrite(ncf,'dr.name',dr.name,[1]);
	addatts(ncf,'name',att);
    nccreate(ncf,'dr.date','Datatype','int32','Dimensions',{'dr.date' length(dr.date)});
	ncwrite(ncf,'dr.date',dr.date,[1]);
	addatts(ncf,'date',att);
    nccreate(ncf,'dr.lat','Datatype','double');
    	ncwrite(ncf,'dr.lat',dr.lat);
    	addatts(ncf,'lat',att);
    nccreate(ncf,'dr.lon','Datatype','double');
    	ncwrite(ncf,'dr.lon',dr.lon);
    	addatts(ncf,'lon',att);
    nccreate(ncf,'dr.ubar','Datatype','double');
    	ncwrite(ncf,'dr.ubar',dr.ubar);
    	addatts(ncf,'ubar',att);
    nccreate(ncf,'dr.vbar','Datatype','double');
    	ncwrite(ncf,'dr.vbar',dr.vbar);
    	addatts(ncf,'vbar',att);

    % zbot dim
    nccreate(ncf,'dr.zbot','Datatype','double','Dimensions',{'zbot' length(dr.zbot)});
	ncwrite(ncf,'dr.zbot',dr.zbot,[1]);
    	addatts(ncf,'zbot',att);
    nccreate(ncf,'dr.ubot','Datatype','double','Dimensions',{'zbot'});
	ncwrite(ncf,'dr.ubot',dr.ubot,[1]);
    	addatts(ncf,'ubot',att);
    nccreate(ncf,'dr.vbot','Datatype','double','Dimensions',{'zbot'});
	ncwrite(ncf,'dr.vbot',dr.vbot,[1]);
    	addatts(ncf,'vbot',att);
    nccreate(ncf,'dr.uerrbot','Datatype','double','Dimensions',{'zbot'});
	ncwrite(ncf,'dr.uerrbot',dr.uerrbot,[1]);
    	addatts(ncf,'uerrbot',att);

    % z_sadcp dim
    nccreate(ncf,'dr.z_sadcp','Datatype','double','Dimensions',{'z_sadcp' length(dr.z_sadcp)});
	ncwrite(ncf,'dr.z_sadcp',dr.z_sadcp,[1]);
    	addatts(ncf,'z_sadcp',att);
    nccreate(ncf,'dr.u_sadcp','Datatype','double','Dimensions',{'z_sadcp'});
	ncwrite(ncf,'dr.u_sadcp',dr.u_sadcp,[1]);
    	addatts(ncf,'u_sadcp',att);
    nccreate(ncf,'dr.v_sadcp','Datatype','double','Dimensions',{'z_sadcp'});
	ncwrite(ncf,'dr.v_sadcp',dr.v_sadcp,[1]);
    	addatts(ncf,'v_sadcp',att);
    nccreate(ncf,'dr.uerr_sadcp','Datatype','double','Dimensions',{'z_sadcp'});
	ncwrite(ncf,'dr.uerr_sadcp',dr.uerr_sadcp,[1]);
    	addatts(ncf,'uerr_sadcp',att);

    % z dim
    nccreate(ncf,'dr.z','Datatype','double','Dimensions',{'z' length(dr.z)});
	ncwrite(ncf,'dr.z',dr.z,[1]);
    	addatts(ncf,'z',att);
    nccreate(ncf,'dr.u','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.u',dr.u,[1]);
    	addatts(ncf,'u',att);
    nccreate(ncf,'dr.v','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.v',dr.v,[1]);
    	addatts(ncf,'v',att);
    nccreate(ncf,'dr.nvel','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.nvel',dr.nvel,[1]);
    	addatts(ncf,'nvel',att);
    nccreate(ncf,'dr.uerr','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.uerr',dr.uerr,[1]);
    	addatts(ncf,'uerr',att);
    nccreate(ncf,'dr.range','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.range',dr.range,[1]);
    	addatts(ncf,'range',att);
    nccreate(ncf,'dr.range_do','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.range_do',dr.range_do,[1]);
    	addatts(ncf,'range_do',att);
    nccreate(ncf,'dr.range_up','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.range_up',dr.range_up,[1]);
    	addatts(ncf,'range_up',att);
    nccreate(ncf,'dr.ts','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.ts',dr.ts,[1]);
    	addatts(ncf,'ts',att);
    nccreate(ncf,'dr.ts_out','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.ts_out',dr.ts_out,[1]);
    	addatts(ncf,'ts_out',att);
    nccreate(ncf,'dr.p','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.p',dr.p,[1]);
    	addatts(ncf,'p',att);
    nccreate(ncf,'dr.ctd_t','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.ctd_t',dr.ctd_t,[1]);
    	addatts(ncf,'ctd_t',att);
    nccreate(ncf,'dr.ctd_s','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.ctd_s',dr.ctd_s,[1]);
    	addatts(ncf,'ctd_s',att);
    nccreate(ncf,'dr.u_do','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.u_do',dr.u_do,[1]);
    	addatts(ncf,'u_do',att);
    nccreate(ncf,'dr.v_do','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.v_do',dr.v_do,[1]);
    	addatts(ncf,'v_do',att);
    nccreate(ncf,'dr.u_up','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.u_up',dr.u_up,[1]);
    	addatts(ncf,'u_up',att);
    nccreate(ncf,'dr.v_up','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.v_up',dr.v_up,[1]);
    	addatts(ncf,'v_up',att);
    nccreate(ncf,'dr.ensemble_vel_err','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.ensemble_vel_err',dr.ensemble_vel_err,[1]);
    	addatts(ncf,'ensemble_vel_err',att);
    nccreate(ncf,'dr.u_shear_method','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.u_shear_method',dr.u_shear_method,[1]);
    	addatts(ncf,'u_shear_method',att);
    nccreate(ncf,'dr.v_shear_method','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.v_shear_method',dr.v_shear_method,[1]);
    	addatts(ncf,'v_shear_method',att);
    nccreate(ncf,'dr.w_shear_method','Datatype','double','Dimensions',{'z'});
	ncwrite(ncf,'dr.w_shear_method',dr.w_shear_method,[1]);
    if existf(dr,'ctd_ss')
	nccreate(ncf,'dr.ctd_ss','Datatype','double','Dimensions',{'z'});
	    ncwrite(ncf,'dr.ctd_ss',dr.ctd_ss,[1]);
	    addatts(ncf,'ctd_ss',att);
    end
    if existf(dr,'ctd_N2')
	nccreate(ncf,'dr.ctd_N2','Datatype','double','Dimensions',{'z'});
	    ncwrite(ncf,'dr.ctd_N2',dr.ctd_N2,[1]);
	    addatts(ncf,'ctd_N2',att);
    end

    % tim dim
    nccreate(ncf,'dr.tim','Datatype','double','Dimensions',{'tim' length(dr.tim)});
	ncwrite(ncf,'dr.tim',dr.tim,[1]);
    	addatts(ncf,'tim',att);
    nccreate(ncf,'dr.tim_hour','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.tim_hour',dr.tim_hour,[1]);
    	addatts(ncf,'tim_hour',att);
    nccreate(ncf,'dr.shiplon','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.shiplon',dr.shiplon,[1]);
    	addatts(ncf,'shiplon',att);
    nccreate(ncf,'dr.shiplat','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.shiplat',dr.shiplat,[1]);
    	addatts(ncf,'shiplat',att);
    nccreate(ncf,'dr.xship','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.xship',dr.xship,[1]);
    	addatts(ncf,'xship',att);
    nccreate(ncf,'dr.yship','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.yship',dr.yship,[1]);
    	addatts(ncf,'yship',att);
    nccreate(ncf,'dr.uship','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.uship',dr.uship,[1]);
    	addatts(ncf,'uship',att);
    nccreate(ncf,'dr.vship','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.vship',dr.vship,[1]);
    	addatts(ncf,'vship',att);
    nccreate(ncf,'dr.zctd','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.zctd',dr.zctd,[1]);
    	addatts(ncf,'zctd',att);
    nccreate(ncf,'dr.wctd','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.wctd',dr.wctd,[1]);
    	addatts(ncf,'wctd',att);
    nccreate(ncf,'dr.uctd','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.uctd',dr.uctd,[1]);
    	addatts(ncf,'uctd',att);
    nccreate(ncf,'dr.vctd','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.vctd',dr.vctd,[1]);
    	addatts(ncf,'vctd',att);
    nccreate(ncf,'dr.xctd','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.xctd',dr.xctd,[1]);
    	addatts(ncf,'xctd',att);
    nccreate(ncf,'dr.yctd','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.yctd',dr.yctd,[1]);
    	addatts(ncf,'yctd',att);
    nccreate(ncf,'dr.uctderr','Datatype','double','Dimensions',{'tim'});
	ncwrite(ncf,'dr.uctderr',dr.uctderr,[1]);
    	addatts(ncf,'uctderr',att);
end % function

%----------------------------------------------------------------------

function [] = add_struct(ncf,snm,struct)

    fname = fieldnames(struct);
    for i=1:length(fname)
	fns = char(fname(i));
	vnm = sprintf('%s.%s',snm,fns);
    
	f = getfield(struct,fns);
	if isstruct(f)
	    error(sprintf('ladcp2cdf:add_struct(%s) substructures are not allowed',vnm));
	end
    
	if ischar(f),		type = 'char';					% define data type
	elseif isnumeric(f),	type = 'double';
	elseif islogical(f),	type = 'int8'; if f, f=1; else, f=0; end;	% logical -> int
	else, error(sprintf('ladcp2cdf:create_var(%s) unsupported type',vnm));
	end

	[nr nc] = size(f);							% define variable
	if (nr*nc < 2)								% scalar
	    nccreate(ncf,vnm,'Datatype',type);
	elseif (nr==1 || nc==1) 						% vector
	    nccreate(ncf,vnm,'Datatype',type,'Dimensions',{vnm nr*nc});
	else									% matrix
	    nccreate(ncf,vnm,'Datatype',type,'Dimensions',{sprintf('%s_c',vnm) nc sprintf('%s_r',vnm) nr});
	end

	if (nr*nc == 1) 							% write data: scalar
	    %disp(sprintf('Writing one %s value to %s...',type,vnm));
	    ncwrite(ncf,vnm,f);
	elseif (nr==1 || nc==1) 						% vector
	    %disp(sprintf('Writing %d %s values to %s...',nr*nc,type,vnm));
	    ncwrite(ncf,vnm,f,[1]);
	elseif (nr*nc > 1)							% matrix
	    %disp(sprintf('Writing %dx%d %s values to %s...',nr,nc,type,vnm));
	    ncwrite(ncf,vnm,f',[1 1]);
        end % if
    end % for
end % function



