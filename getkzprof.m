function  [dk,pk]=getkzprof(dr,pk,N2);
% function  [dk,pk]=getkzprof(dr,pk,N2);
%
% get profile of Kz
% input dr: result from LADCP inversion
%       pk: control parameter
%       N2: stratification  = -1/rho0 dsig/dz [1/s^2]
%

%======================================================================
%                    G E T K Z P R O F . M 
%                    doc: Thu Jun 24 01:10:55 2004
%                    dlm: Wed Jan  7 16:27:30 2009
%                    (c) 2004 ladcp@
%                    uE-Info: 20 1 NIL 0 0 72 0 2 8 NIL ofnI
%======================================================================

% CHANGES BY ANT:
%  Jun 24, 2004: - typo in plot legend (cast TOO shallow)
%  Jan  7, 2009: - tightened use of exist()

if ~exist('pk','var')
 pk.top='default';
end

% get vertical resolution
dz=mean(diff(dr.z));

% desired vertical range over which to compute spectrum
pk=setdefv(pk,'zrange',500);

% number of grid points to compute spectra
imax=ceil(pk.zrange/dz);
imax=2.^round(log2(imax));
pk.kz_imax=imax;

disp([' use dz:',int2str(dz*pk.kz_imax),' m long segments'])
dk.zrange=dz*pk.kz_imax;

% number of grid points to slide down
pk=setdefv(pk,'kz_ioff',5);

% use also shearbased profile if available
pk=setdefv(pk,'use_shear',1);

% first depth index
iz=(1:pk.kz_imax)+pk.kz_ioff;

% loop over profile
m=0;
while max(iz)<length(dr.z)
   m=m+1;
   pk.iuse=iz;
   % get kz for this segment
   if nargin>2
    dkk=getkz(dr,pk,N2);
   else
    dkk=getkz(dr,pk);
   end
   % save results to matrix
   dk.Kz(m)=dkk.K_z;
   dk.z(m)=mean(dr.z(iz));
   dk.N(m)=dkk.N;
   dk.Eps(m)=dkk.Eps;
   dk.Kz_err(m)=dkk.K_z_err;
   dk.Kz_min(m)=dkk.K_z_min;
   dk.Kz_max(m)=dkk.K_z_max;
   iz=iz+pk.kz_ioff;
   dk.N2_type=dkk.N2_type;
end

dk.name=dr.name;
dk.date=dr.date;
dk.lat=dr.lat;
dk.lon=dr.lon;

figure(8)
clf,
orient tall
streamer([dr.name,' Figure 8']);
if existf(dk,'z')==0, 

title([' CAST TOO SHALLOW ']) 

else

semilogx(dk.Kz,-dk.z,'-')
hold on
semilogx(dk.Kz_min,-dk.z,'--')
semilogx(dk.Kz_max,-dk.z,'--')
ax=axis;
ylabel('depth [m]')
xlabel('K_z [m^2 s^{-1}]')
title([' Experimental K_z profile   ',dk.N2_type]) 

end
