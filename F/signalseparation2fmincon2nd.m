
function curvedistance = signalseparation2fmincon2nd(ng1, ng2, examplecurve, coeffs2)

% coeffs(coeffs<0)=0;
% coeffs(coeffs>1)=1;

G1_coeff = coeffs2(1);
G2_coeff = coeffs2(2);

examplecurve = examplecurve / max(examplecurve(1));
examplecurve = max(min(examplecurve,1),0);

combined_curve = G1_coeff * ng1 + G2_coeff * ng2;
curve_difference = combined_curve - examplecurve;

delta = 1;
curvedistance = mean(delta*curve_difference.^2 + (1-delta) * abs(curve_difference)) ...
    + 0.1 * abs(coeffs2(1))...%加入了强稀疏约束，强制某一通道信号变少
    - 0 * abs(coeffs2(1)-coeffs2(2))%加入了对两个通道差值的约束，强迫他们变大
% lambda = 0.5;
% curvedistance = mean(curve_difference.^2) + lambda * norm(coeffs,2)^2;%加入正则化
% curvedistance = mean(curve_difference.^2);
% curvedistance =  mean(curve_difference);