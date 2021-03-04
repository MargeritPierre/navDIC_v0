%% LCURVE ASSOCIATED TO A TIKHNOV REGULARIZATION
% TEST FOR A NUMBER OF VALUES OF BETA
BETA = 10.^linspace(3,7,10) ;

% Constant objects
H0 = weight*Hess(validDOF,validDOF)*weight ;
Hr = CONS ;
J0 = weight*dr_da(validDOF) ;
Jr = - CONS*DU ;

nTest = numel(BETA) ;
tic
R0 = NaN(nTest,1) ;
Rr = NaN(nTest,1) ;
N = NaN(nTest,1) ;
for bb = 1:nTest
    A = H0+BETA(bb)*Hr ;
    b = J0+BETA(bb)*Jr ;
    x = A\b ;
    N(bb) = norm(x) ;
    R0(bb) = norm(H0*x-J0) ;
    Rr(bb) = norm(Hr*x-Jr) ;
end
toc

clf ;
plot3(R0/max(R0),Rr/max(Rr),BETA)
set(gca,'xscale','log','yscale','log','zscale','log')
xlabel('$\| H_0 \cdot x - j_0 \|$','interpreter','latex')
ylabel('$\| H_r \cdot x - j_r \|$','interpreter','latex')
zlabel('$\beta$','interpreter','latex')



