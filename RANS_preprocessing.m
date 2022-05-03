clear
close all

case_name = 'r150_ccl_np12';
C = colororder;
blk = read_grid(case_name);
NB = length(blk);
figure(1)
hold on
for i=1:NB
    plot(blk{i}.x,blk{i}.y,'color',C(mod(i-1,7)+1,:));
    plot(blk{i}.x',blk{i}.y','color',C(mod(i-1,7)+1,:));
end
axis equal

figure(2)
hold on
for i=1:NB
    x = blk{i}.x;
    y = blk{i}.y;
    ni = size(x,1);
    nj = size(x,2);
    plot(x(ni,:),y(ni,:),'k')
    plot(x(1,:)',y(1,:)','k')
    plot(x(:,1),y(:,1),'k')
    plot(x(:,nj),y(:,nj),'k')
end
axis equal
figure(3)
hold on
skip = 4;
for i=1:NB
    x = blk{i}.x;
    y = blk{i}.y;
    ni = size(x,1);
    nj = size(x,2);
    for j=1:skip:ni
        plot(x(j,:),y(j,:),'color',C(mod(i-1,7)+1,:))
    end
    for j=1:skip:nj
        plot(x(:,j),y(:,j),'color',C(mod(i-1,7)+1,:))
    end
    plot(x(end,:),y(end,:),'color',C(mod(i-1,7)+1,:))
    plot(x(:,end),y(:,end),'color',C(mod(i-1,7)+1,:))
end
set(gca,'FontSize',14)
axis equal
%%
figure(4)
hold on
skip = 4;
for i=1:NB
    x = blk{i}.x;
    y = blk{i}.y;
    ni = size(x,1);
    nj = size(x,2);
    for j=1:skip:ni
        if (j==1) || (j==ni)
            plot(x(j,:),y(j,:),'r','LineWidth',1)
        else
            plot(x(j,:),y(j,:),'k')
        end
    end
    for j=1:skip:nj
        if (j==1) || (j==ni)
            plot(x(j,:),y(j,:),'r','LineWidth',1)
        else
            plot(x(:,j),y(:,j),'k')
        end
    end
    plot(x(1,:),y(1,:),'r','LineWidth',1)
    plot(x(end,:),y(end,:),'r','LineWidth',1)
    plot(x(:,1),y(:,1),'r','LineWidth',1)
    plot(x(:,end),y(:,end),'r','LineWidth',1)
end

x = blk{6}.x(:,1);
y = blk{6}.y(:,1);
(x(end)-x(1))/length(x)

set(gca,'FontSize',14)
axis equal

base_folder = cd;
temp_slash = '/'; if ispc, temp_slash = '\'; end
%write_plot3d_2d(blk,[base_folder temp_slash case_name temp_slash case_name '.x']);
