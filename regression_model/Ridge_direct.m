function w_estimate = Ridege_direct(Output, Dic, lambda)
%%
% This is a VANILA implementation of the follwing paper
%
% The algorithm solves the inverse problem
%   $y = Aw + \xi$
%
% ============= Inputs ===============
% y                                    : output, in the paper $y(t) =
%                                           (x(t+\delta)-x(t))/\delta$;
% A                                    : dictionary matrix;
% lambda                         : the tradeoff parameter you should use,
%                                          basically, it is proportional to
%                                          the invese of the variance, e.g. 1;
% MAXITER                    : maximum number of the iterations, e.g. 5

% ============= Reference =============
% W. Pan, Y. Yuan, J. Goncalves, and G.-B. Stan,
% A Sparse Bayesian Approach to the Iden- tification of Nonlinear State-Space Systems,
% IEEE Transaction on Automatic Control, 2015 (to appear). arXiv:1408.3549
% http://arxiv.org/abs/1408.3549
%
% ============= Author =============
%  Wei Pan (w.pan11@imperial.ac.uk, panweihit@gmail.com)
%
% ============= Version =============
%   1.0 (Sep ?, 2012)
%%

delta = 1e-4;

[M,N]=size(Dic);

W = (lambda * eye(N) + Dic'*Dic )\Dic'* Output ;

W(W./norm(W)<sqrt(delta))=0;
w_estimate = W ;




