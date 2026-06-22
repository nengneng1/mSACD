function curvedistance = signalseparation2fminconBG(ng1, ng2, examplecurve, coeffs)

% coeffs(coeffs<0)=0;
% coeffs(coeffs>1)=1;

G1_coeff = coeffs(1);
G2_coeff = coeffs(2);
bg = coeffs(3);

examplecurve = examplecurve / examplecurve(1);
examplecurve = max(min(examplecurve,1),0);

combined_curve = bg + G1_coeff * ng1 + G2_coeff * ng2;
curve_difference = abs(combined_curve - examplecurve);
curvedistance = sqrt(mean(curve_difference.^2));
% curvedistance = mean(curve_difference);