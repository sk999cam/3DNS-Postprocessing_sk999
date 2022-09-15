function write_input_files(casename,blk,next_block,next_patch,corner,bcs,gas,solver,topology,nkproc)
nargin;
if nargin < 10 || ~exist('nkproc','var')
    nkproc = ceil(solver.nk/solver.npp);
    fprintf('k procs: %d\n', nkproc)
end

if ~isfield(gas,'gamma')
    gas.gamma = gas.gam;
end
npp = solver.npp;
fprintf('Using npp = %d\n', npp)
NB = length(blk.x);
ncorner = length(corner);   

dir = fullfile(pwd,casename);
fprintf('Writing to directory: %s\n',dir)

if(~exist(dir,'dir'))
mkdir(dir);
end

% Now write header for new input file
% GPU input
fidin = fopen(fullfile(dir,'input_gpu.txt'),'w');
fprintf(fidin,'%d %d\n', [NB,1]);
   
nprocs = 0;

for ib=1:NB
    
x = blk.x{ib};
y = blk.y{ib};
[ni,nj]=size(x);   

im_next_block = next_block{ib}.im;
ip_next_block = next_block{ib}.ip;
jm_next_block = next_block{ib}.jm;
jp_next_block = next_block{ib}.jp;

im_next_patch = next_patch{ib}.im;
ip_next_patch = next_patch{ib}.ip;
jm_next_patch = next_patch{ib}.jm;
jp_next_patch = next_patch{ib}.jp;
    

%nkproc = 1;
%niproc = 1;%max([1,floor(ni/(4*npp))]);
%njproc = 1;%max([1,floor(nj/(4*npp))]);
%nprocs = nprocs + niproc*njproc*nkproc;

fprintf(fidin,'%d %d %d\n', [ni nj solver.nk]);
  
%fprintf(fidin,'%d %d %d\n', [niproc njproc nkproc]);    
   
      
   im = (im_next_block==0)*im_next_patch;
   ip = (ip_next_block==0)*ip_next_patch;
   
   jm = (jm_next_block==0)*jm_next_patch;
   jp = (jp_next_block==0)*jp_next_patch;
   
   
   fprintf(fidin,'%d %d %nid %d\n', [im ip jm jp]);
   
   if(im==0)
   fprintf(fidin,'%d %d\n', [im_next_block im_next_patch]);
   end
   
   if(ip==0)
   fprintf(fidin,'%d %d\n', [ip_next_block ip_next_patch]);
   end
   
   if(jm==0)
   fprintf(fidin,'%d %d\n', [jm_next_block jm_next_patch]);
   end
   
   if(jp==0)
   fprintf(fidin,'%d %d\n', [jp_next_block jp_next_patch]);
   end
   
   NI(ib) = ni;
   NJ(ib) = nj;
   
end

% write corners
fprintf(fidin,'%d\n', ncorner);
for n=1:ncorner
    cor_type(n) = 0;
end

for n=1:ncorner
fprintf(fidin,'%d %d\n', [corner{n}.Nb, cor_type(n)]);
for nb=1:corner{n}.Nb
ib = corner{n}.block{nb};
ic = corner{n}.i{nb};
jc = corner{n}.j{nb};
if(ic>1); ic = NI(ib); end
if(jc>1); jc = NJ(ib); end
fprintf(fidin,'%d %d %d\n', [ib ic jc]);
end

end

fprintf(fidin,'%d\n',1); % 1 block group
fprintf(fidin,'%d\n',NB); % NB blocks in group
for ib=1:NB
fprintf(fidin,'%d ',ib); % blocks in group
end

% 
% write rest of file
% check these are ok for your case!
    
    %nsteps nwrite ncut 
    fprintf(fidin,'\n%d %d %d\n', [solver.niter solver.nwrite solver.ncut]);
    
    % cfl, filter coefficient
    fprintf(fidin,'%f %f\n', [solver.cfl solver.sigma]);
    
     % Toin poin pext vref alpha_in pitch_in aturb (not used) ilength (not used) g_z
    fprintf(fidin,'%f %f %f %f %f %f %f %d %d %f\n', [bcs.Toin bcs.Poin bcs.pexit bcs.vin bcs.alpha bcs.gamma bcs.aturb bcs.ilength bcs.radprof bcs.g_z]);
    
    % gamma cp mu_ref sutherlands constants prandtl no.
    fprintf(fidin,'%f %f %12.5e %f %f %f\n', [gas.gamma gas.cp gas.mu_ref gas.mu_tref gas.mu_cref gas.pr]);
    
    % span, spanwise grid expansion factor
    fprintf(fidin,'%f %f\n', [solver.span solver.fexpan]);
    
    % restart, statistics
    fprintf(fidin,'%d %d\n', [solver.irestart solver.istats]);

fclose(fidin);



% Now write header for new input file
% CPU input
fidin = fopen(fullfile(dir,'input_cpu.txt'),'w');
fprintf(fidin,'%d\n', [NB]);
   
nprocs = 0;

for ib=1:NB
    
x = blk.x{ib};
y = blk.y{ib};
[ni,nj]=size(x);
nk = solver.nk;

im_next_block = next_block{ib}.im;
ip_next_block = next_block{ib}.ip;
jm_next_block = next_block{ib}.jm;
jp_next_block = next_block{ib}.jp;

im_next_patch = next_patch{ib}.im;
ip_next_patch = next_patch{ib}.ip;
jm_next_patch = next_patch{ib}.jm;
jp_next_patch = next_patch{ib}.jp;
    


niproc = ceil(ni/npp);
njproc = ceil(nj/npp);
nprocs = nprocs + niproc*njproc*nkproc;

fprintf(fidin,'%d %d %d\n', [ni nj nk]);
  
fprintf(fidin,'%d %d %d\n', [niproc njproc nkproc]);    
   
      
   im = (im_next_block==0)*im_next_patch;
   ip = (ip_next_block==0)*ip_next_patch;
   
   jm = (jm_next_block==0)*jm_next_patch;
   jp = (jp_next_block==0)*jp_next_patch;
   
   
   fprintf(fidin,'%d %d %d %d\n', [im ip jm jp]);
   
   if(im==0)
   fprintf(fidin,'%d %d\n', [im_next_block im_next_patch]);
   end
   
   if(ip==0)
   fprintf(fidin,'%d %d\n', [ip_next_block ip_next_patch]);
   end
   
   if(jm==0)
   fprintf(fidin,'%d %d\n', [jm_next_block jm_next_patch]);
   end
   
   if(jp==0)
   fprintf(fidin,'%d %d\n', [jp_next_block jp_next_patch]);
   end
   
   NI(ib) = ni;
   NJ(ib) = nj;
   
end

% write corners
fprintf(fidin,'%d\n', ncorner);
for n=1:ncorner
    cor_type(n) = 0;
end

for n=1:ncorner
fprintf(fidin,'%d %d\n', [corner{n}.Nb, cor_type(n)]);
for nb=1:corner{n}.Nb
ib = corner{n}.block{nb};
ic = corner{n}.i{nb};
jc = corner{n}.j{nb};
if(ic>1); ic = NI(ib); end
if(jc>1); jc = NJ(ib); end
fprintf(fidin,'%d %d %d\n', [ib ic jc]);
end

end
% 
% write rest of file
% check these are ok for your case!
    
    %nsteps nwrite ncut 
    fprintf(fidin,'%d %d %d\n', [solver.niter solver.nwrite solver.ncut]);
    
    % cfl, filter coefficient, ifsplit, ifsat, if_LES, if
    fprintf(fidin,'%f %f %d %d %d %d %d\n', [solver.cfl, solver.sigma solver.ifsplit solver.ifsat solver.ifLES 0 0]);
    
     % Toin poin pext vref alpha_in pitch_in aturb (not used) ilength (not used) g_z
    fprintf(fidin,'%f %f %f %f %f %f %f %f %d %d\n', [bcs.Toin bcs.Poin bcs.pexit bcs.vin bcs.alpha bcs.cax bcs.aturb bcs.lturb bcs.ilength bcs.radprof]);
    
    % gamma cp mu_ref sutherlands constants prandtl no.
    fprintf(fidin,'%f %f %12.5e %f %f %f\n', [gas.gamma gas.cp gas.mu_ref gas.mu_tref gas.mu_cref gas.pr]);
    
    % span, spanwise grid expansion factor
    fprintf(fidin,'%f %f\n', [solver.span solver.fexpan]);
    
    % restart, statistics
    fprintf(fidin,'%d %d\n', [solver.irestart solver.istats]);
    
    % inlet groups
    fprintf(fidin,'%d\n', [1]);

    if topology == 1
        fprintf(fidin,'%d\n', [2]);   
        fprintf(fidin,'%d\n%d\n', [1 2]);
    elseif topology == 2
        fprintf(fidin,'%d\n', [3]);   
        fprintf(fidin,'%d\n%d\n%d\n', [1 2 3]);
    end
     
    % stability stuff   
    fprintf(fidin,'%d %d %d %d %d %f %f\n', [solver.istability 150 177 183 2 5.0 10000.0]);
   
    % inlet boundary layer stuff and adiabatic walls
    fprintf(fidin,'%d %f %f', [0 1e5 bcs.twall]);


fclose(fidin);
fprintf('Total processors: %d\n', nprocs)

end