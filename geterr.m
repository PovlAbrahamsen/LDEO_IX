%======================================================================
%                    G E T E R R . M 
%                    doc: Wed Jun 30 23:24:51 2004
%                    dlm: Fri Mar  5 15:48:05 2010
%                    (c) 2004 ladcp@
%                    uE-Info: 91 17 NIL 0 0 72 2 2 8 NIL ofnI
%======================================================================

% MODIFICATIONS BY ANT:
%  Jun 30, 2004: - BUG: bin numbering was wrong for asymmetric up/down
%	                bin setup
%  Jul  5, 2004: - added comments to debug depth mapping
%  Jul 16, 2004: - added global variable skip_figure_3 to workaround
%		   linux matlab bug
%  Oct  7, 2008: - extensively modified procfig 3 for version IX_6

function l=geterr(dr,d,iplot)
% function l=geterr(dr,d,iplot)
% returns predicitons of U_ocean and
% U_ctd on the raw data grid
% 
% CTD velocity
if nargin<3, iplot=1; end

tim=dr.tim;
tim(1)=-1e30;
tim(end)=1e30;

uctd=-interp1q(tim',dr.uctd',d.time_jul');
vctd=-interp1q(tim',dr.vctd',d.time_jul');

[ib,it]=size(d.ru);

wm=medianan(d.rw,3);
wz=gradient(-d.z,d.time_jul*24*3600);
l.ru_ctd=meshgrid(uctd,1:ib)+d.weight*0;
l.rv_ctd=meshgrid(vctd,1:ib)+d.weight*0;
l.rw_ctd=meshgrid(wm,1:ib)+d.weight*0;
l.rw_ctd_z=meshgrid(wz,1:ib)+d.weight*0;
if existf(d,'wctd')
 l.rw_ctd_p=meshgrid(d.wctd,1:ib)+d.weight*0;
end

% OCEAN velocity

z=-d.izm+d.ru*0;
dz=diff(d.izm(:,1))';

ii=find(z>=min(dr.z) & z<=max(dr.z));

uoce=interp1q(dr.z,dr.u,z(ii));
voce=interp1q(dr.z,dr.v,z(ii));

[prof,bin]=meshgrid(1:it,1:ib);

l.ru_oce=full(sparse(bin(ii),prof(ii),uoce));
l.rv_oce=full(sparse(bin(ii),prof(ii),voce));
l.ru_oce(ib,it)=NaN;
l.rv_oce(ib,it)=NaN;
l.ru_oce=l.ru_oce+d.weight*0;
l.rv_oce=l.rv_oce+d.weight*0;
ii=find(l.ru_oce==0 & l.rv_oce==0);
l.ru_oce(ii)=NaN;
l.rv_oce(ii)=NaN;

% ocean velocity as a function of depth and time

					% ib is number of bins
					% it is number of times (super ensembles)
itm=meshgrid(1:it,1:ib);		% each of ib rows of itm contains 1:it

					% d.izm contains for each time (colums),
					% list of absolute depths for each bin
dzdo=mean(abs(diff(d.izm(d.izd,1))));	% dzdo contains sound-speed corrected
					% mean bin length of downlooker at surface
					% NB: at depth, bins are smaller, because
					%     of increased soundspeed!

if length(d.izu)>1			% uplooker bin length
 dzup=mean(abs(diff(d.izm(d.izu,1))));
else
 dzup=dzdo;
end
dz=min([dzdo dzup]);			% dz is min bin length near surface
iz=-(d.izm/dz);				% iz is d.izm with depth coordinate given
					% as distance from surface, measured in 
					% near-surface bin lengths 

					% d.ru contains super-ensemble velocities
					% dr.z contains output depth grid
ii=find(isfinite(d.ru) & iz>0 & iz<max(dr.z)/dz);
					% ii contains indices (valid for d.ru,
					% d.izm, iz, ...) with valid velocities,
					% inside the output depth grid

ij=find( iz>0 & iz<max(dr.z)/dz);	% ij contains same as ii but also for
					% invalid velocities

if abs(dzup-dzdo)>dzup*0.1
 disp([' sorry dz not constant loop ',int2str(length(ii)),' elements'])
 for j=1:length(ii)
  iiz=ceil(iz(ii(j)));
  iit=itm(ii(j));
  l.u_oce(iiz,iit)=d.ru(ii(j))-l.ru_ctd(ii(j));
  l.v_oce(iiz,iit)=d.rv(ii(j))-l.rv_ctd(ii(j));
  l.w_oce(iiz,iit)=d.rw(ii(j))-l.rw_ctd(ii(j));
  l.w_oce_z(iiz,iit)=d.rw(ii(j))-l.rw_ctd_z(ii(j));
  if existf(l,'rw_ctd_p')
   l.w_oce_p(iiz,iit)=d.rw(ii(j))-l.rw_ctd_p(ii(j));
  end
  if existf(d,'tg')
   l.tg_oce(iiz,iit)=d.tg(ii(j));
  end

  l.u_ocean(iiz,iit)=l.ru_oce(ii(j));
  l.v_ocean(iiz,iit)=l.rv_oce(ii(j));

  l.u_adcp(iiz,iit)=d.ru(ii(j));
  l.v_adcp(iiz,iit)=d.rv(ii(j));
 end
else % uplooker and downlooker bin sizes are equal
 l.u_oce=full(sparse(ceil(iz(ii)),itm(ii),d.ru(ii)-l.ru_ctd(ii)));
 l.v_oce=full(sparse(ceil(iz(ii)),itm(ii),d.rv(ii)-l.rv_ctd(ii)));
 l.w_oce=full(sparse(ceil(iz(ii)),itm(ii),d.rw(ii)-l.rw_ctd(ii)));
 l.w_oce_z=full(sparse(ceil(iz(ii)),itm(ii),d.rw(ii)-l.rw_ctd_z(ii)));
 if existf(l,'rw_ctd_p')
  l.w_oce_p=full(sparse(ceil(iz(ii)),itm(ii),d.rw(ii)-l.rw_ctd_p(ii)));
 end
 if existf(d,'tg')
  l.tg_oce=full(sparse(ceil(iz(ij)),itm(ij),d.tg(ij)));
 end

 l.u_ocean=full(sparse(ceil(iz(ii)),itm(ii),l.ru_oce(ii)));
 l.v_ocean=full(sparse(ceil(iz(ii)),itm(ii),l.rv_oce(ii)));

 l.u_adcp=full(sparse(ceil(iz(ii)),itm(ii),d.ru(ii)));
 l.v_adcp=full(sparse(ceil(iz(ii)),itm(ii),d.rv(ii)));
end

ik=find(l.u_oce==0 & l.v_oce==0);
l.u_oce(ik)=NaN;
l.v_oce(ik)=NaN;
l.w_oce(ik)=NaN;
l.w_oce_z(ik)=NaN;
if existf(l,'rw_ctd_p')
 l.w_oce_p(ik)=NaN;
end
l.u_adcp(ik)=NaN;
l.v_adcp(ik)=NaN;
if existf(d,'tg')
 ik=find(l.tg_oce==0);
 l.tg_oce(ik)=NaN;
end

[lz,lt]=size(l.u_oce);
l.itv=1:lt;

l.z_oce=([1:lz]-.5)*dz;
l.u_oce_m=meannan(l.u_oce');
l.v_oce_m=meannan(l.v_oce');

l.u_oce_s=stdnan(l.u_oce');
l.v_oce_s=stdnan(l.v_oce');

l.ru_err=d.ru-l.ru_oce-l.ru_ctd;
l.rv_err=d.rv-l.rv_oce-l.rv_ctd;

l.izm=d.izm;

[lz,lt]=size(l.ru_err);
l.itv2=1:lt;

if iplot

% blank out shallow/deep estimates
ii=find(iz<0 | iz>max(dr.z)/dz);
d.ru(ii)=nan;
d.rv(ii)=nan;

global skip_figure_3;
if isempty(skip_figure_3)

   figure(3)
   clf
   orient landscape
   
   colormap(polarmap(21));
   
   subplot(231)
   ib=1:size(l.ru_err,1);
   ib=ib-length(d.izu);
   tmp = l.ru_err; tmp(isnan(tmp)) = 0;
   pcolorn(l.itv2,-ib,tmp), shading flat
   fac=meannan(l.u_oce_s);
   fac=max([fac,1e-2]);
   caxis([-3 3]*fac)
   colorbar
   xlabel('Super Ensemble #');
   ylabel('Bin #');
   title(sprintf('U-err std: %.03f',meannan(stdnan(l.ru_err'))))
   
   subplot(232)
   plot(meannan(l.ru_err')',-ib)
   set(gca,'XLim',[-0.05 0.05]);
   set(gca,'Ylim',[-ib(end) -ib(1)]);
   set(gca,'Xtick',[-0.04:0.02:0.04]);
   grid
   xlabel('Residual [m/s]');
   ylabel('Bin #');
   title('mean(U-err)')
   
   subplot(233)
   tmp = l.u_oce; tmp(isnan(tmp)) = 0;
   pcolorn(l.itv,-l.z_oce,tmp), shading flat
   ca = caxis;
   if abs(ca(1)) > abs(ca(2))
    caxis([-abs(ca(1)) abs(ca(1))]);
   else 
    caxis([-abs(ca(2)) abs(ca(2))]);
   end
   if existf(dr,'zbot')
    hold on
    plot(-d.z+d.hbot,'.k')
    ax=axis;
    ax(4)=maxnan([-d.z+d.hbot,ax(4)]);
    axis(ax);
   end
   colorbar
   xlabel('Ensemble #');
   ylabel('Depth [m]');
   title('U_{oce}')
   
   subplot(234)
   tmp = l.rv_err; tmp(isnan(tmp)) = 0;
   pcolorn(l.itv2,-ib,tmp), shading flat
   fac=meannan(l.v_oce_s);
   fac=max([fac,1e-2]);
   caxis([-3 3]*fac)
   colorbar
   xlabel('Super Ensemble #');
   ylabel('Bin #');
   title(sprintf('V-err std: %.03f',meannan(stdnan(l.rv_err'))))
   
   subplot(235)
   plot(meannan(l.rv_err')',-ib)
   set(gca,'XLim',[-0.05 0.05]);
   set(gca,'Ylim',[-ib(end) -ib(1)]);
   set(gca,'Xtick',[-0.04:0.02:0.04]);
   grid
   xlabel('Residual [m/s]');
   ylabel('Bin #');
   title('mean(V-err)')

   subplot(236)
   tmp = l.v_oce; tmp(isnan(tmp)) = 0;
   pcolorn(l.itv,-l.z_oce,tmp), shading flat
   ca = caxis;
   if abs(ca(1)) > abs(ca(2))
    caxis([-abs(ca(1)) abs(ca(1))]);
   else 
    caxis([-abs(ca(2)) abs(ca(2))]);
   end
   if existf(dr,'zbot')
    hold on
    plot(-d.z+d.hbot,'.k')
    ax=axis;
    ax(4)=maxnan([-d.z+d.hbot,ax(4)]);
    axis(ax);
   end
   colorbar
   xlabel('Ensemble #');
   ylabel('Depth [m]');
   title('U_{oce}')
   
   streamer([dr.name,'  Figure 3']);
end % of ~skip_figure_3

% reset colormap
figure(11)
colormap(jet(128))

end


%======================================================================
%                    P O L A R M A P . M 
%                    doc: Tue Oct  7 11:03:28 2008
%                    dlm: Tue Oct  7 11:13:04 2008
%                    (c) 2008 A.M. Thurnherr
%                    uE-Info: 21 51 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

function map = polarmap(n)

if nargin<1, n=129, end;	% same as for jet()

map = ones(n,3);

firstred  = ceil(n/2) + 1;
lastblue = floor(n/2);

map([1:lastblue],1) = [0:lastblue-1]'/lastblue;
map([1:lastblue],2) = [0:lastblue-1]'/lastblue;
map([firstred:end],2) = [lastblue-1:-1:0]'/lastblue;
map([firstred:end],3) = [lastblue-1:-1:0]'/lastblue;
