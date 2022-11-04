function write_turbid_input(turbid, casename)
    
	lsx = turbid.lsx;
	lsy = turbid.lsy;
	lsz = turbid.lsz;
	nfx = turbid.nfx;
	nfy = turbid.nfy;
	nfz = turbid.nfz;
    ni  =turbid.ni;
    nj  =turbid.nj;
    nk  =turbid.nk;

	basedir = pwd;
	path = fullfile(basedir,casename);
    if ~exist(path,'dir')
		mkdir(path);
    end

	turbid_path = fullfile(path,'turbid_files');
	if ~exist(turbid_path,'dir')
		mkdir(turbid_path);
    end
    

	turb_file = 'turbid.txt'  ;
	file_path = fullfile(turbid_path,turb_file);
	fprintf('Writing turbid input file for case %s\n', casename);
		
	% Now write input file for turbid

	f = fopen(file_path, 'w');
	fprintf(f,'%d %d %d\n', [ni,nj,nk]);
	fprintf(f,'%d %d %d\n', [lsx,lsy,lsz]);
	fprintf(f,'%d %d %d\n', [nfx,nfy,nfz]);
	fclose(f)
end
    
    