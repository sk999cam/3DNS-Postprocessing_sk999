clear
close all

case_name = 'r150_pitch_7';
C = colororder;
blk = read_grid(case_name);
ni_inlet = 15;
len_inlet = 1;

dx = blk{1}.x(2,1)-blk{1}.x(1,1);
x0 = blk{1}.x(1,1);
x_inlet=x0*ones(1,ni_inlet);
r0 = 1;
r1 = 10;
r = 0.5*(r0+r1);
while (r-r0)/r > 0.0001
    S = dx*(1-r^(ni_inlet-1))/(1-r);
    if S > len_inlet
        r1 = r;
    else
        r0 = r;
    end
    r = (r0+r1)/2;
end
for i=1:ni_inlet-2
    x_inlet(ni_inlet-i) = x_inlet(ni_inlet+1-i) - dx*r^(i-1);
end
x_inlet(1) = x0-len_inlet;
y_inlet = blk{1}.y(1,:);
[X,Y] = meshgrid(x_inlet,y_inlet);
X = X';
Y = Y';
blk{length(blk)+1}.x = X;
blk{length(blk)}.y = Y;
y_inlet = blk{2}.y(1,:);
[X,Y] = meshgrid(x_inlet,y_inlet);
X = X';
Y = Y';
blk{length(blk)+1}.x = X;
blk{length(blk)}.y = Y;

NB = length(blk);

%%
figure(1)
hold on
for i=1:NB
    plot(blk{i}.x,blk{i}.y,'color',C(mod(i-1,7)+1,:));
    plot(blk{i}.x',blk{i}.y','color',C(mod(i-1,7)+1,:));
end
axis equal

base_folder = cd;
temp_slash = '/'; if ispc, temp_slash = '\'; end
%write_plot3d_2d(blk,[base_folder temp_slash case_name temp_slash case_name 'long_inlet.x']);
