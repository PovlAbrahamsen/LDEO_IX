function [dk,pk]=getkz(dr,pk,N2,dens)
% function [dk,pk]=getkz(dr,pk,N2,dens)
%
% compute vertical diffusivity from velocity profile
% 
% follow Polzin et al. JAOTEC 2002 page 205
%
% first compute velocity shear spectra
%  also compute displacement spectra
% 
%  M. Visbeck Aug 2002
%  modified Feb 2004

if nargin<2, pk.top='default'; end

% variables used
%  N0 = reference Brunt Vaisala Frequency,   1/s
%   b = reference depth scale for Garret and Munk spectrum
%   j = parameter for Garrent and Munk spectrum
%

% set some default values
% default N0, j and b for Garret and Munk spectrum
pk=setdefv(pk,'N0',5.24e-3);
pk=setdefv(pk,'j',3);
pk=setdefv(pk,'b',1300);

if nargin<3, 
 % if N2 data do not exist use GM N^2 profile
 %     N(z) = N0 * e ^(-(z-z0)/b)
 N=pk.N0*exp(-(abs(dr.z)-200)/pk.b); 
 ii=find(abs(dr.z)<1500);
 N(ii)=N(ii(end));
 %N=pk.N0+dr.z*0;
 N2=N.^2;
 dk.N2_type='use GM reference';
else
 dk.N2_type='from CTD data';
end

if nargin<4, 
 dens=N2*nan;
 dk.dens_type='no density data';
end

% if a reference N is given make sure it has the correct length
if length(N2)~=length(dr.z)
 N2=mean(N2)+0*dr.z;
end

% compute N ignore negative values
N=sqrt(abs(N2));

% Use complex velocity profiles
u=dr.u+i*dr.v;
p=dr.z;
pk=setdefv(pk,'use_shear',0);
if pk.use_shear & existf(dr,'u_shear_method')
 u_s=dr.u_shear_method+i*dr.v_shear_method;
else
 pk.use_shear=0;
end

% find typical dz
dz=median(diff(dr.z));

% find maximal power of 2
i2max=fix(log2(length(dr.z)));
pk=setdefv(pk,'i2max',i2max);

% find maximal length of profile that can be used,
% limit to middle of profile
iuse=fix((length(u)-2^pk.i2max)/2)+[1:2.^pk.i2max];
pk=setdefv(pk,'iuse',iuse);

% remove trend or polynome from data
pk=setdefv(pk,'polyfit',1);
pk.polyfit=max(pk.polyfit,0);
c=polyfit(p(pk.iuse),u(pk.iuse),pk.polyfit);
u=u(pk.iuse)-polyval(c,p(pk.iuse));

if pk.use_shear 
  c=polyfit(p(pk.iuse),u_s(pk.iuse),pk.polyfit);
  u_s=u_s(pk.iuse)-polyval(c,p(pk.iuse));
end

N=N(pk.iuse);
dens=dens(pk.iuse);

% reference stratification
N2=median(N).^2;

pk=setdefv(pk,'i2fft',pk.i2max);


% Compute spectrum of observed shear data
nn=2.^pk.i2fft;
% SPECTRUM
fs=2*pi/dz;                % sampling frequ
df=fs/nn;               % frequ. interval
nn2=fix(nn/2);
kv=(1:(nn2-1)).*df; % nn/2 frequ
n0=nn/2;

% u and v spectrum
us=spectrum(real(u).*N,nn,n0);
vs=spectrum(imag(u).*N,nn,n0);
if length(isfinite(dens))==length(dens)
 hs=spectrum(dens,nn,n0);
end
 
%us=sqrt(us.^2+vs.^2);
us=us+vs;
Us=us(2:nn2)*dz.*kv.^2;
if pk.use_shear
 u_ss=spectrum(real(u_s).*N,nn,n0);
 v_ss=spectrum(imag(u_s).*N,nn,n0);
% u_ss=sqrt(u_ss.^2+v_ss.^2);
 u_ss=u_ss+v_ss;
 Uss=u_ss(2:nn2)*dz.*kv.^2;
 Us=mean([Us;Uss]);
 dk.Us_shear=Uss;
 dk.Us_inverse=Us;
end
Us=Us./mean(N).^2;
% correction due to Kees Veth (have to check at some point)
Us=Us.*fs/2; 


% make spectral correction for LADCP sampling
% Kunze et al. (Jaotec 2001)
pk=setdefv(pk,'dzr',dz/(2*pi));
pk=setdefv(pk,'sinc_exp',12);

dk.T_cor=sinc(kv*pk.dzr).^pk.sinc_exp;

% correct shear spectrum for LADCP sampling error
Us=Us./dk.T_cor;

% Compute GM shear spectrum
[UGM,UsGM,HGM]=velgm(kv,sqrt(N2),pk.N0,pk.j,pk.b);

% GET epsilon
% set shear/strain ratio to 1
% this could be done using the displacement spectra
% not yet implemented
pk=setdefv(pk,'fRw',1);
% set epsilon 0
pk=setdefv(pk,'e0',7.8e-10);

efac=pk.e0.*N2./pk.N0.^2*pk.fRw;
e=efac*Us.^2./UsGM.^2;

dkv=mean(diff(kv));
Usint=cumsum(Us)/N2*dkv;

% GET K_z
% set gamma to 0.2  gamma=Rf/(1-Rf)
pk=setdefv(pk,'ga',0.2);

% K_z = gamma e / N^2
Kz=pk.ga.*e./N2;
% Kz if shear was at GM level
Kz0=pk.ga.*efac./N2;

% find best Kz within reasonable range
pk=setdefv(pk,'kvav',[0.002 0.09]);
pk=setdefv(pk,'median',1);
ii=find(kv>=pk.kvav(1) & kv<=pk.kvav(2));

% limit range to intergal < 0.7
i2=sum(Usint(ii)<=0.7);

% dont go below 3 points
i2=max(i2,3);

if length(ii)>2
 ii=ii(1:i2);
 if pk.median
  Eps=median(e(ii));
  Kza=median(Kz(ii));
 else
  Eps=mean(e(ii));
  Kza=mean(Kz(ii));
 end
 Kza(2)=min(Kz(ii));
 Kza(3)=max(Kz(ii));
 Kza(4)=std(Kz(ii))/sqrt(length(ii));
 dk.kv_av=[kv(ii(1)) kv(ii(end))];
else
 Kza=[NaN NaN NaN];
 dk.kv_av=pk.kvav;
 Eps=NaN;
end

dk.UsGM=UsGM;
dk.Us=Us;
dk.kv=kv;
dk.Eps=Eps;
dk.K_z=Kza(1);
dk.K_z_min=Kza(2);
dk.K_z_max=Kza(3);
dk.K_z_err=Kza(4);
dk.N=sqrt(N2);

figure(8)
clf,
loglog(kv,UsGM/N2,'-r')
hold on 
loglog(kv,Us./N2,'-b','linewidth',2)
loglog(kv,(UsGM./N2)*sqrt(Kza(1)./Kz0),'--b')
ax=axis;
loglog(dk.kv_av(1)*[1,1],ax(3:4),'-k')
loglog(dk.kv_av(2)*[1,1],ax(3:4),'-k')
ylabel('S[V_z / N](k_z) 1/(rad/m)')
xlabel('k_z (rad/m)')
title([' Mean Kz: ',num2str(Kza(1)*1e5),' +/- ',num2str(Kza(4)*1e5),' 10^{-5}  m^2/s  Epsilon: ',...
       num2str(Eps),' W/kg'])
legend('GM shear','LADCP shear','fit to LADCP')


% ----------------------------------------
function [Uk,Usk,Hk,kv0]=velgm(kv,N,N0,j,b)
%
% from Gregg and Kunze (JGR 1991 page 16.717)
% see also Cairns and Williams (JGR 1976 page 1943) 
%
% INPUT:
%  kv: vertical wave number rad/m
%  N: stratification in rad/s  0.00524 = 3 cycles/hour
%  N0: ref stratification in rad/s  0.00524 = 3 cycles/hour
%  j: mode number (3)
%  b: thermocline depth (1300 m)
%
% OUTPUT:
%  Uk: horizontal velocity spectra [m^2 s^-2/(rad m^-1)]
%  Usk: horizontal shear spectra   [s^-2 /(rad m^-1)]
%  Hk: vertical displacement spectra [m^2/(rad m^-1)]
%  kv0: reference wave number      [rad m^-1]
%

if nargin<5, b=1300; end
if nargin<4, j=3; end
if nargin<3, N0=0.00524; end
if nargin<2, N=N0; end

% vertial reference wave number
kv0=pi*j/b*(N./N0);

% Energy level of GM spectrum (no dimension)
E=6.3e-5;

% Horizontal Velocity Spectrum 
Uk=3*E*b.^3*N0.^2/(2*j*pi)./(1+kv./kv0).^2;

% Horizontal Velocity Shear Spectrum 
Usk=kv.^2 .* Uk;

% Vertical Displacement Spectrum 
Hk = E*b.^3*(N0./N).^2/(2*j*pi)./(1+kv./kv0).^2;
