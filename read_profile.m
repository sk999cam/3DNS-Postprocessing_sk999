function [xprof,yprof,pitch,stag]=read_profile(file_name, inorm)

fid = fopen(file_name,'r');  
N=fscanf(fid, '%d',1);
pitch=fscanf(fid, '%f',1);
for i=1:N
A=fscanf(fid, '%f %f',2);
xprof(i) = A(1);
yprof(i) = A(2);
end
fclose(fid);

[xLE,iLE] = min(xprof);
yLE = yprof(iLE);

[xTE,iTE] = max(xprof);
yTE = yprof(iTE);

flip_prof = sign(yprof(iLE)-yprof(mod(iLE+1,N)));

cax = xTE-xLE;
stag = atan((yTE-yLE)/(xTE-xLE))*180/pi;

if flip_prof == 1
    xprof = flip(xprof);
    yprof = flip(yprof);
end

if inorm
    yprof = (yprof-yLE)/cax;
    xprof = (xprof-xLE)/cax;
    pitch = pitch/cax;
end

return