function [Rd,RS,RFv11,RFv12,RFv21,RFv22,RFv33] = residual_MS(mp,dt,T,F,Fe,Fv,q,S,Fn,Fvn,qn,Sn)
etaD0 = mp.etaD0;
ks = mp.ks;
etaS = mp.etaS;
etaN = mp.etaN;
a0 = mp.a0;
a1 = mp.a1;
a2 = mp.a2;
Tni = mp.Tni;
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
Fv11 = Fv(1,1);
Fv12 = Fv(1,2);
Fv21 = Fv(2,1);
Fv22 = Fv(2,2);
Fv33 = Fv(3,3);
Fn11 = Fn(1,1);
Fn12 = Fn(1,2);
Fn21 = Fn(2,1);
Fn22 = Fn(2,2);
Fn33 = Fn(3,3);
Fvn11 = Fvn(1,1);
Fvn12 = Fvn(1,2);
Fvn21 = Fvn(2,1);
Fvn22 = Fvn(2,2);
Fvn33 = Fvn(3,3);

Rd = q - qn + dt*(cos(q)^2*((F12 - Fn12)/(2*F22) - (F12*(F11 - Fn11))/(2*F11*F22)) + sin(q)^2*((F12 - Fn12)/(2*F22) - (F12*(F11 - Fn11))/(2*F11*F22))) - (3*S*dt*exp((2*((muneq*(((Fe11*Fe21 + Fe12*Fe22)*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + ((Fe11*Fe21 + Fe12*Fe22)*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*S*cos(q)*sin(q)*(Fe11^2 + Fe12^2))/((2*S + 1)*(S - 1)) + (3*S*cos(q)*sin(q)*(Fe21^2 + Fe22^2))/((2*S + 1)*(S - 1))))/2 + (mueq*((F12*F22*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (F12*F22*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*F22^2*S*cos(q)*sin(q))/((2*S + 1)*(S - 1)) + (3*S*cos(q)*sin(q)*(F11^2 + F12^2))/((2*S + 1)*(S - 1))))/2)^2 + ((muneq*((2*(Fe11^2 + Fe12^2)*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (6*S*cos(q)*sin(q)*(Fe11*Fe21 + Fe12*Fe22))/((2*S + 1)*(S - 1)) - 2))/2 + (muneq*((2*Fe33^2)/(S - 1) + 2))/2 + (mueq*((2*(F11^2 + F12^2)*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (6*F12*F22*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))))/2 + mueq/(F11^2*F22^2*(S - 1)))^2 + ((muneq*((2*(Fe21^2 + Fe22^2)*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (6*S*cos(q)*sin(q)*(Fe11*Fe21 + Fe12*Fe22))/((2*S + 1)*(S - 1)) - 2))/2 + (muneq*((2*Fe33^2)/(S - 1) + 2))/2 + (mueq*((2*((3*S*sin(q)^2)/(2*S + 1) - 1)*F22^2)/(S - 1) + (6*F12*S*cos(q)*sin(q)*F22)/((2*S + 1)*(S - 1))))/2 + mueq/(F11^2*F22^2*(S - 1)))^2)^(1/2)/ks)*(cos(q)*(sin(q)*(mueq*(F11^2 + F12^2) + muneq*(Fe11^2 + Fe12^2)) - cos(q)*(muneq*(Fe11*Fe21 + Fe12*Fe22) + F12*F22*mueq)) + sin(q)*(sin(q)*(muneq*(Fe11*Fe21 + Fe12*Fe22) + F12*F22*mueq) - cos(q)*(mueq*F22^2 + muneq*(Fe21^2 + Fe22^2)))))/(etaD0*(2*S + 1)*(S - 1));
RS = S - Sn + (dt*(S^3*a2 - S^2*a1 + S*a0*(T - Tni) + (mueq*(F11^2 - sin(q)*((F22^2*sin(q)*(6*S^2 + 3))/(2*S + 1)^2 + (F12*F22*cos(q)*(6*S^2 + 3))/(2*S + 1)^2) - cos(q)*((cos(q)*(6*S^2 + 3)*(F11^2 + F12^2))/(2*S + 1)^2 + (F12*F22*sin(q)*(6*S^2 + 3))/(2*S + 1)^2) + F12^2 + F22^2 + 1/(F11^2*F22^2)))/(2*(S - 1)^2) + (muneq*(Fe11^2 + Fe12^2 + Fe21^2 + Fe22^2 + Fe33^2 - cos(q)*((sin(q)*(Fe11*Fe21 + Fe12*Fe22)*(6*S^2 + 3))/(2*S + 1)^2 + (cos(q)*(6*S^2 + 3)*(Fe11^2 + Fe12^2))/(2*S + 1)^2) - sin(q)*((cos(q)*(Fe11*Fe21 + Fe12*Fe22)*(6*S^2 + 3))/(2*S + 1)^2 + (sin(q)*(6*S^2 + 3)*(Fe21^2 + Fe22^2))/(2*S + 1)^2)))/(2*(S - 1)^2)))/etaS;
RFv11 = Fv11 - Fvn11 - (Fv11*dt*muneq*(Fe11*((Fe11*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe21*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) + Fe21*((Fe21*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe11*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) - 1))/etaN - (Fv21*dt*muneq*(Fe12*((Fe11*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe21*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) + Fe22*((Fe21*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe11*S*cos(q)*sin(q))/((2*S + 1)*(S - 1)))))/etaN;
RFv12 = Fv12 - Fvn12 - (Fv12*dt*muneq*(Fe11*((Fe11*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe21*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) + Fe21*((Fe21*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe11*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) - 1))/etaN - (Fv22*dt*muneq*(Fe12*((Fe11*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe21*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) + Fe22*((Fe21*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe11*S*cos(q)*sin(q))/((2*S + 1)*(S - 1)))))/etaN;
RFv21 = Fv21 - Fvn21 - (Fv11*dt*muneq*(Fe11*((Fe12*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe22*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) + Fe21*((Fe22*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe12*S*cos(q)*sin(q))/((2*S + 1)*(S - 1)))))/etaN - (Fv21*dt*muneq*(Fe12*((Fe12*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe22*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) + Fe22*((Fe22*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe12*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) - 1))/etaN;
RFv22 = Fv22 - Fvn22 - (Fv12*dt*muneq*(Fe11*((Fe12*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe22*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) + Fe21*((Fe22*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe12*S*cos(q)*sin(q))/((2*S + 1)*(S - 1)))))/etaN - (Fv22*dt*muneq*(Fe12*((Fe12*((3*S*cos(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe22*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) + Fe22*((Fe22*((3*S*sin(q)^2)/(2*S + 1) - 1))/(S - 1) + (3*Fe12*S*cos(q)*sin(q))/((2*S + 1)*(S - 1))) - 1))/etaN;
RFv33 = Fv33 - Fvn33 + (Fv33*dt*muneq*(Fe33^2/(S - 1) + 1))/etaN;

% W12 = dt*(cos(q)^2*((F12 - Fn12)/(2*F22) - (F12*(F11 - Fn11))/(2*F11*F22)) + sin(q)^2*((F12 - Fn12)/(2*F22) - (F12*(F11 - Fn11))/(2*F11*F22)))
end