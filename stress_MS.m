function sigma = stress_MS(mp,F,Fe,q,S)
mueq = mp.mueq;
muneq = mp.muneq;

F11 = F(1,1);
F12 = F(1,2);
F22 = F(2,2);
Fe11 = Fe(1,1);
Fe12 = Fe(1,2);
Fe21 = Fe(2,1);
Fe22 = Fe(2,2);
Fe33 = Fe(3,3);
 
s11 = (muneq*((2*(Fe11^2 + Fe12^2)*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (6*S*cos(q)*sin(q)*(Fe11*Fe21 + Fe12*Fe22))/((2*S + 1)*(S - 1)) - 2))/2 + (muneq*((2*Fe33^2)/(S - 1) + 2))/2 + (mueq*((2*(F11^2 + F12^2)*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (6*F12*F22*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))))/2 + mueq/(F11^2*F22^2*(S - 1));
s12 = (muneq*(((Fe11*Fe21 + Fe12*Fe22)*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + ((Fe11*Fe21 + Fe12*Fe22)*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*S*cos(q)*sin(q)*(Fe11^2 + Fe12^2))/((2*S + 1)*(S - 1)) + (3*S*cos(q)*sin(q)*(Fe21^2 + Fe22^2))/((2*S + 1)*(S - 1))))/2 + (mueq*((F12*F22*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (F12*F22*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*F22^2*S*cos(q)*sin(q))/((2*S + 1)*(S - 1)) + (3*S*cos(q)*sin(q)*(F11^2 + F12^2))/((2*S + 1)*(S - 1))))/2;
s22 = (muneq*((2*(Fe21^2 + Fe22^2)*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (6*S*cos(q)*sin(q)*(Fe11*Fe21 + Fe12*Fe22))/((2*S + 1)*(S - 1)) - 2))/2 + (muneq*((2*Fe33^2)/(S - 1) + 2))/2 + (mueq*((2*((3*S*sin(q)^2)/(2*S + 1) - 1)*F22^2)/(S - 1) + (6*F12*S*cos(q)*sin(q)*F22)/((2*S + 1)*(S - 1))))/2 + mueq/(F11^2*F22^2*(S - 1));

sigma = zeros(3,3);
sigma(1,1) = s11;
sigma(1,2) = s12;
sigma(2,1) = s12;
sigma(2,2) = s22;
end