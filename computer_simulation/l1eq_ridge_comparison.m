function l1eq_ridge_comparison(input, pool_range, tur, lambda_l2, repeat, outdir, core)
%on server

%{
    input:
	   input: directory and the name of the input table, header = TRUE, row = sample, column = gene, for example: input = '/home/singlecell/dataset1.txt'
           pool_range: the range of pool, e.g. pool_range = [10:40] or [15,20,25,30,35,40,45]
           tur: turbulance, e.g. tur = 0.2
           lambda_l2: the tradeoff parameter for Ridge regression-based method,basically, it is proportional to the invese of the variance, e.g. 10^(-2) 
           repeat: the repeat times
           outdir: the directory of output
           core: core used during calculation
		   
	output: 
           Several data matrix include:
	       A datamat, a cormedianmat, a cormeanmat
%}

X = importdata(input);
X = X.data;
cell = size(X,1)
datamat = zeros(8,size(pool_range,2));
corcellmedianmat = zeros(2*(repeat+1),size(pool_range,2));
corcellmeanmat = zeros(2*(repeat+1),size(pool_range,2));
expression = zeros(2*(repeat+1),size(pool_range,2));

col = 0;
row = 1;

for pool = pool_range
    
    col = col+1;
    corcellmedianmat(1,col) = pool;
    corcellmeanmat(1,col) = pool;
    expression(1,col) = pool;
    corcellmedianmat(repeat+2,col) = pool;
    corcellmeanmat(repeat+2,col) = pool;
    expression(repeat+2,col) = pool;
    
end 

p = parpool(core);
p.IdleTimeout = 100000000;

    for t = 1:repeat
      
    row = row+1;
    col = 0;
    
    for pool = pool_range
    
      col = col+1;
      phi = randi([0,1],pool,cell);
      e = unifrnd(1-tur,1+tur,pool,cell);
      ephi = e.*phi;     
      
      result_l1eq = l1eq_simulation_full(X,ephi,phi);
      corcellmedianmat(row,col) = result_l1eq(3,1);
      corcellmeanmat(row,col) = result_l1eq(4,1);
      expression(row,col) = result_l1eq(5,1);
      screen = strcat('l1_eq simulation finished! pool = ',num2str(pool),', repeat = ',num2str(t),', time = ', datestr(datetime('now')));
      disp(screen)

      result_l2 = ridge_simulation_full(X, ephi, phi, lambda_l2);
      corcellmedianmat(row+(repeat+1),col) = result_l2(3,1);
      corcellmeanmat(row+(repeat+1),col) = result_l2(4,1);
      expression(row+(repeat+1),col) = result_l2(5,1);
      screen = strcat('Ridge simulation finished! pool = ',num2str(pool),', repeat = ',num2str(t),', time = ', datestr(datetime('now')));
      disp(screen)
      
screen = strcat('pool = ',num2str(pool),', repeat = ',num2str(t),' finished! time = ', datestr(datetime('now')));
disp(screen)

    filename_corcellmedian = strcat(outdir,'simulation_tur',num2str(tur),'_repeat',num2str(repeat),'_corcellmedianmat.mat');
    save(filename_corcellmedian,'corcellmedianmat');
    filename_corcellmean = strcat(outdir,'simulation_tur',num2str(tur),'_repeat',num2str(repeat),'_corcellmeanmat.mat');
    save(filename_corcellmean,'corcellmeanmat');
    filename_express = strcat(outdir,'simulation_tur',num2str(tur),'_repeat',num2str(repeat),'_expression.mat');
    save(filename_express,'expression');
    

    end
    

    end

col = 0;
for pool = pool_range
    col = col+1;
    datamat(1,col) = pool;
    datamat(5,col) = pool;
    
    datamat(2,col) = nanmedian(corcellmedianmat(2:(repeat+1),col)) ;
    datamat(3,col) = nanmean(corcellmeanmat(2:(repeat+1),col)) ;
    datamat(4,col) = nanmean(expression(2:(repeat+1),col));
    
    
    datamat(6,col) = nanmedian(corcellmedianmat((repeat+3):2*(repeat+1),col)) ;
    datamat(7,col) = nanmean(corcellmeanmat((repeat+3):2*(repeat+1),col)) ;
    datamat(8,col) = nanmean(expression((repeat+3):2*(repeat+1),col)) ;
    
end 

filename = strcat(outdir,'simulation_tur',num2str(tur),'_repeat',num2str(repeat),'_lambdal2_',num2str(lambda_l2),'_datamat.mat');
save(filename,'datamat');

delete(gcp); 

