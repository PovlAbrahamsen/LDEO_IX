% README for LADCP software
%
% Martin Visbeck
disp(' LADCP software ')
disp(' ')
disp(' M. Visbeck  LDEO  ')
disp(' http://www.ldeo.columbia.edu/~visbeck/ladcp ')
disp(' ')
a=which('readme');
ii=find(a=='/');
ladcpdir='';
if length(ii)>0
 ladcpdir=a(1:max(ii));
end
disp([' You m-files are installed in directory :',ladcpdir])
disp(' ')
disp(' Much of the parameter are documented in DEFAULT.M ')
disp(' press RETURN to display it ')

pause
more on
eval(' type default ')
pause
clc

disp(' The data are processed by  LAPROC.M ')
disp(' press RETURN to display it ')

pause
eval(' type laproc ')

pause
clc

disp(' Copy any of the load?.m files to your local directory') 
eval(['ls ',ladcpdir,'loaddata/load*.m'])

pause
clc

disp(' Edit any of the demo?.m files to your needs') 
eval(['ls ',ladcpdir,'process/d*.m'])

pause
clc

disp('% parameter structure arrays used')
disp(' ')
disp('p.*   main parameter')
disp('ps.*  inverse solution parameter')
disp('f.* file names ')
disp(' ')
disp('% data structure arrays used')
disp(' ')
disp('d.*   raw data')
disp('di.*  super ensemble data')
disp('dr.*  results from inversion')
disp('de.*  matrix data used by inversion')
disp('ds.*  results from shear based method')
disp('da.*  description of profile for NETCDF file save')
more off

