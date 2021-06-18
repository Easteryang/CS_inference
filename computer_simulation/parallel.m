function l1eq_ridge_comparison(input, blockSize, ratio, outdir, core)
%on server

%{
    input:
	    input: directory and the name of the input table, header = TRUE, row = sample, column = gene, for example: input = '/home/singlecell/big_dataset.txt'
	    blockSize: the number of cells per block, for example: 500
	    ratio: the number of pools per block / the number of cells per block, belongs to the range of (0,1), for example: 0.5
	    outdir: the directory of output
        core: core used during calculation
		   
	output: 
           Several data matrices include:
	       A recoverdatamat, a correlationdatamat
%}

cellData=importdata(input);
addpath(genpath('/path that contains ridge_simulation.m'));

blockNum = floor(size(cellData,2)/blockSize);
poolSize = floor(blockSize*ratio);

Rarray = {};
for i = 1:blockNum
  curData = cellData( : , ((i-1)*blockSize+1):(i*blockSize) );
  pool = poolSize;
  Rarray{i} = ridge_simulation(curData, pool, core);
  % disp(i);
  % disp('Done!');
end

remainSize = size(cellData,2)-blockNum*blockSize;
if ( remainSize > 0 )
	i = blockNum+1;
	curData = cellData( : , (blockNum*blockSize+1):size(cellData,2) );
	poolSize = floor(remainSize*ratio);
	pool = poolSize
	Rarray{i} = ridge_simulation(curData, pool, core);
	% disp(i);
	% disp('Done!');
end

recoverData = [];
for j = 1:length(Rarray)
    recoverData = [recoverData;Rarray{j}];
end

recoverData = recoverData';

filename = strcat(outdir,'block',num2str(blockSize),'_ratio',num2str(ratio),'_recoverData.mat');
save(filename,'recoverData','-v7.3');

m = corr(cellData,recoverData);

filename = strcat(outdir,'block',num2str(blockSize),'_ratio',num2str(ratio),'_correlation.mat');
save(filename,'m','-v7.3');