function l1eq_experiment(input,inputphi,outdir, core)
%on server

%{
    input:
	   input: directory and the name of the input table, header = FALSE, row = gene, column = pool, for example: input = '/home/singlecell/experiment.txt'
           inputphi: normalized measurement matrix, row = pool, column = sample, for example: inputphi = '/home/singlecell/phi_normalization.mat'
           outdir: the directory of output
           core: core used during calculation
		   
	output: 
           Several data matrix include:
	       A recover datamat, row = gene, column = sample
%}

addpath(genpath('/path/regression_model'))

X = importdata(input);
X = X';
pool= size(X,1)
gene = size(X,2)

%Calling multiple cores
p = parpool(core);
p.IdleTimeout = 100000000;

phi = importdata(inputphi);
cell = size(phi,2)

recoverX = [];

%Multiple threads running simultaneously
parfor i = 1:size(X,2)
    y = X(:,i); 
    x0 = phi'*y
    xp = l1eq_pd(x0, phi, [], y);
    recoverX(:,i) = xp;
end

recoverX(recoverX<0) = 0;

filename_recoverX = strcat(outdir,'basis_simulation','_recoverx.mat');
save(filename_recoverX,'recoverX');

delete(gcp); 
