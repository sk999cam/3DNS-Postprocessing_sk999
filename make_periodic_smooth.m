function [blk] = make_periodic_smooth(blk,next_block,next_patch,nbnow,pitch)

ng = 1;

up.xi = 0;
up.yi = 0;
dn.xi = 0;
dn.yi = 0;

up.xj = 0;
up.yj = 0;
dn.xj = 0;
dn.yj = 0;

ib = nbnow;

im_next_block = next_block{ib}.im;
ip_next_block = next_block{ib}.ip;
jm_next_block = next_block{ib}.jm;
jp_next_block = next_block{ib}.jp;

im_next_patch = next_patch{ib}.im;
ip_next_patch = next_patch{ib}.ip;
jm_next_patch = next_patch{ib}.jm;
jp_next_patch = next_patch{ib}.jp;



if(im_next_block~=0)

if(im_next_patch==1)
yoffset = 0.0;
if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end  

xup = blk{im_next_block}.x(2,:);
yup = blk{im_next_block}.y(2,:)+yoffset;

blk{ib}.x(1,:) = (xup + blk{ib}.x(2,:))*0.5;
blk{ib}.y(1,:) = (yup + blk{ib}.y(2,:))*0.5;
blk{im_next_block}.x(1,:) = blk{ib}.x(1,:);
blk{im_next_block}.y(1,:) = blk{ib}.y(1,:)-yoffset;
end

if(im_next_patch==2)
yoffset = 0.0;
if( (blk{ib}.y(1,1) - blk{im_next_block}.y(end,1)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,1) - blk{im_next_block}.y(end,1)) < -0.1*pitch ); yoffset =-pitch; end          
xup = blk{im_next_block}.x(end-1,:);
yup = blk{im_next_block}.y(end-1,:)+yoffset;

blk{ib}.x(1,:) = (xup + blk{ib}.x(2,:))*0.5;
blk{ib}.y(1,:) = (yup + blk{ib}.y(2,:))*0.5;
blk{im_next_block}.x(end,:) = blk{ib}.x(1,:);
blk{im_next_block}.y(end,:) = blk{ib}.y(1,:)-yoffset;
end

if(im_next_patch==3)
yoffset = 0.0;
if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,1)') >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,1)') < -0.1*pitch ); yoffset =-pitch; end          
xup = blk{im_next_block}.x(:,2)';
yup = blk{im_next_block}.y(:,2)'+yoffset;


blk{ib}.x(1,:) = (xup + blk{ib}.x(2,:))*0.5;
blk{ib}.y(1,:) = (yup + blk{ib}.y(2,:))*0.5;
blk{im_next_block}.x(:,1) = blk{ib}.x(1,:)';
blk{im_next_block}.y(:,1) = blk{ib}.y(1,:)'-yoffset;
end

if(im_next_patch==4)
yoffset = 0.0;
if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,end)') >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,1) - blk{im_next_block}.y(1,end)') < -0.1*pitch ); yoffset =-pitch; end              
xup = blk{im_next_block}.x(:,end-1)';
yup = blk{im_next_block}.y(:,end-1)'+yoffset;


blk{ib}.x(1,:) = (xup + blk{ib}.x(2,:))*0.5;
blk{ib}.y(1,:) = (yup + blk{ib}.y(2,:))*0.5;
blk{im_next_block}.x(:,end) = blk{ib}.x(1,:)';
blk{im_next_block}.y(:,end) = blk{ib}.y(1,:)'-yoffset;
end

up.xi = xup;
up.yi = yup;

end

if(ip_next_block~=0)

if(ip_next_patch==1)
yoffset = 0.0;
if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end         
xdn = blk{ip_next_block}.x(2,:);
ydn = blk{ip_next_block}.y(2,:)+yoffset;


blk{ib}.x(end,:) = (xdn + blk{ib}.x(end-1,:))*0.5;
blk{ib}.y(end,:) = (ydn + blk{ib}.y(end-1,:))*0.5;
blk{ip_next_block}.x(1,:) = blk{ib}.x(end,:);
blk{ip_next_block}.y(1,:) = blk{ib}.y(end,:)-yoffset;
end

if(ip_next_patch==2)
yoffset = 0.0;
if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(end,1)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(end,1)) < -0.1*pitch ); yoffset =-pitch; end          
xdn = blk{ip_next_block}.x(end-1,:);
ydn = blk{ip_next_block}.y(end-1,:)+yoffset;

blk{ib}.x(end,:) = (xdn + blk{ib}.x(end-1,:))*0.5;
blk{ib}.y(end,:) = (ydn + blk{ib}.y(end-1,:))*0.5;
blk{ip_next_block}.x(end,:) = blk{ib}.x(end,:);
blk{ip_next_block}.y(end,:) = blk{ib}.y(end,:)-yoffset;
end

if(ip_next_patch==3)
yoffset = 0.0;
if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,1)') >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,1)') < -0.1*pitch ); yoffset =-pitch; end              
xdn = blk{ip_next_block}.x(:,2)';
ydn = blk{ip_next_block}.y(:,2)'+yoffset;

blk{ib}.x(end,:) = (xdn + blk{ib}.x(end-1,:))*0.5;
blk{ib}.y(end,:) = (ydn + blk{ib}.y(end-1,:))*0.5;
blk{ip_next_block}.x(:,1) = blk{ib}.x(end,:)';
blk{ip_next_block}.y(:,1) = blk{ib}.y(end,:)'-yoffset;
end

if(ip_next_patch==4)
yoffset = 0.0;    
if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,end)') >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(end,1) - blk{ip_next_block}.y(1,end)') < -0.1*pitch ); yoffset =-pitch; end              
xdn = blk{ip_next_block}.x(:,end-1)';
ydn = blk{ip_next_block}.y(:,end-1)'+yoffset;

blk{ib}.x(end,:) = (xdn + blk{ib}.x(end-1,:))*0.5;
blk{ib}.y(end,:) = (ydn + blk{ib}.y(end-1,:))*0.5;
blk{ip_next_block}.x(:,end) = blk{ib}.x(end,:)';
blk{ip_next_block}.y(:,end) = blk{ib}.y(end,:)'-yoffset;
end

dn.xi = xdn;
dn.yi = ydn;

end




if(jm_next_block~=0)

if(jm_next_patch==1)
yoffset = 0.0;
if( (blk{ib}.y(1,1)' - blk{jm_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,1)' - blk{jm_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end      
xup = blk{jm_next_block}.x(2,:)';
yup = blk{jm_next_block}.y(2,:)'+yoffset;

blk{ib}.x(:,1) = (xup + blk{ib}.x(:,2))*0.5;
blk{ib}.y(:,1) = (yup + blk{ib}.y(:,2))*0.5;
blk{jm_next_block}.x(1,:) = blk{ib}.x(:,1)';
blk{jm_next_block}.y(1,:) = blk{ib}.y(:,1)'-yoffset;
end

if(jm_next_patch==2)
yoffset = 0.0;
if( (blk{ib}.y(1,1)' - blk{jm_next_block}.y(end,1)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,1)' - blk{jm_next_block}.y(end,1)) < -0.1*pitch ); yoffset =-pitch; end          
xup = blk{jm_next_block}.x(end-1,:)';
yup = blk{jm_next_block}.y(end-1,:)'+yoffset;

blk{ib}.x(:,1) = (xup + blk{ib}.x(:,2))*0.5;
blk{ib}.y(:,1) = (yup + blk{ib}.y(:,2))*0.5;
blk{jm_next_block}.x(end,:) = blk{ib}.x(:,1)';
blk{jm_next_block}.y(end,:) = blk{ib}.y(:,1)'-yoffset;
end

if(jm_next_patch==3)
yoffset = 0.0;
if( (blk{ib}.y(1,1) - blk{jm_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,1) - blk{jm_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end   
xup = blk{jm_next_block}.x(:,2);
yup = blk{jm_next_block}.y(:,2)+yoffset;

blk{ib}.x(:,1) = (xup + blk{ib}.x(:,2))*0.5;
blk{ib}.y(:,1) = (yup + blk{ib}.y(:,2))*0.5;
blk{jm_next_block}.x(:,1) = blk{ib}.x(:,1);
blk{jm_next_block}.y(:,1) = blk{ib}.y(:,1)-yoffset;
end

if(jm_next_patch==4)
yoffset = 0.0;
if( (blk{ib}.y(1,1) - blk{jm_next_block}.y(1,end)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,1) - blk{jm_next_block}.y(1,end)) < -0.1*pitch ); yoffset =-pitch; end
xup = blk{jm_next_block}.x(:,end-1);
yup = blk{jm_next_block}.y(:,end-1)+yoffset;

blk{ib}.x(:,1) = (xup + blk{ib}.x(:,2))*0.5;
blk{ib}.y(:,1) = (yup + blk{ib}.y(:,2))*0.5;
blk{jm_next_block}.x(:,end) = blk{ib}.x(:,1);
blk{jm_next_block}.y(:,end) = blk{ib}.y(:,1)-yoffset;
end

up.xj = xup;
up.yj = yup;

end


if(jp_next_block~=0)

if(jp_next_patch==1)
yoffset = 0.0;    
if( (blk{ib}.y(1,end)' - blk{jp_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,end)' - blk{jp_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end          
xdn = blk{jp_next_block}.x(2,:)';
ydn = blk{jp_next_block}.y(2,:)'+yoffset;

blk{ib}.x(:,end) = (xdn + blk{ib}.x(:,end-1))*0.5;
blk{ib}.y(:,end) = (ydn + blk{ib}.y(:,end-1))*0.5;
blk{jp_next_block}.x(1,:) = blk{ib}.x(:,end)';
blk{jp_next_block}.y(1,:) = blk{ib}.y(:,end)'-yoffset;
end

if(jp_next_patch==2)
yoffset = 0.0;
if( (blk{ib}.y(1,end)' - blk{jp_next_block}.y(end,1)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,end)' - blk{jp_next_block}.y(end,1)) < -0.1*pitch ); yoffset =-pitch; end              
xdn = blk{jp_next_block}.x(end-1,:)';
ydn = blk{jp_next_block}.y(end-1,:)'+yoffset;

blk{ib}.x(:,end) = (xdn + blk{ib}.x(:,end-1))*0.5;
blk{ib}.y(:,end) = (ydn + blk{ib}.y(:,end-1))*0.5;
blk{jp_next_block}.x(end,:) = blk{ib}.x(:,end)';
blk{jp_next_block}.y(end,:) = blk{ib}.y(:,end)'-yoffset;
end

if(jp_next_patch==3)
yoffset = 0.0;
if( (blk{ib}.y(1,end) - blk{jp_next_block}.y(1,1)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,end) - blk{jp_next_block}.y(1,1)) < -0.1*pitch ); yoffset =-pitch; end   
xdn = blk{jp_next_block}.x(:,2);
ydn = blk{jp_next_block}.y(:,2)+yoffset;

blk{ib}.x(:,end) = (xdn + blk{ib}.x(:,end-1))*0.5;
blk{ib}.y(:,end) = (ydn + blk{ib}.y(:,end-1))*0.5;
blk{jp_next_block}.x(:,1) = blk{ib}.x(:,end);
blk{jp_next_block}.y(:,1) = blk{ib}.y(:,end)-yoffset;
end

if(jp_next_patch==4)
yoffset = 0.0;
if( (blk{ib}.y(1,end) - blk{jp_next_block}.y(1,end)) >  0.1*pitch ); yoffset = pitch; end
if( (blk{ib}.y(1,end) - blk{jp_next_block}.y(1,end)) < -0.1*pitch ); yoffset =-pitch; end
xdn = blk{jp_next_block}.x(:,end-1);
ydn = blk{jp_next_block}.y(:,end-1)+yoffset;

blk{ib}.x(:,end) = (xdn + blk{ib}.x(:,end-1))*0.5;
blk{ib}.y(:,end) = (ydn + blk{ib}.y(:,end-1))*0.5;
blk{jp_next_block}.x(:,end) = blk{ib}.x(:,end);
blk{jp_next_block}.y(:,end) = blk{ib}.y(:,end)-yoffset;
end

dn.xj = xdn;
dn.yj = ydn;

end






return