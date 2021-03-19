function resultmat = ridge_simulation_full(X, ephi, phi, lambda_l2)
%The full version adds three new data, the median and average of the correlation coefficients by gene, and the signal-to-interference ratio (SIR)

addpath(genpath('/path/regression_model'))

resultmat = zeros(8,1);
recoverX = zeros(size(X));

parfor i = 1:size(X,2)
    y = ephi*X(:,i);          
	recover = Ridge_direct(y,phi,lambda_l2);
    recoverX(:,i) = recover;
end

recoverX(recoverX<0) = 0;
idxNan = isnan(recoverX);
idxN = (sum(idxNan)>0);
recoverX1 = recoverX(:,~idxN);
X1 = X(:,~idxN);

cor = corr(recoverX1',X1')';

resultmat(1,1) = norm(recoverX1-X1,'fro'); 
resultmat(2,1) = norm(recoverX1-X1,'fro')/norm(X1,'fro');
resultmat(3,1) = nanmedian(diag(cor)); 
resultmat(4,1) = nanmean(diag(cor));  

express_recover = sum(recoverX1>0.001);
express1 = sum(express_recover>0);
resultmat(5,1) = express1;

cor2 = zeros(1,size(X1,2));
SIR = zeros(1,size(X1,2));
for i = 1:size(X1,2)
    cor2(1,i)= corr(recoverX1(:,i),X1(:,i));
    SIR(1,i) = 20*log10( norm(X1(:,i)) / norm(X1(:,i)-recoverX1(:,i)) );
end
resultmat(6,1) = nanmedian(cor2);
resultmat(7,1) = nanmean(cor2);
resultmat(8,1) = nanmean(SIR);
