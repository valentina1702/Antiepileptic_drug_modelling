% EN.580.640 - Final Project
% Sophia Nehs, Valentina Dsouza, Caroline Ghio, Christianne Chua, Shruthi Bare
% Levetiracetam Equations

function dydt = LEV_eqns(t,y,p)
%
% Equations describing levetiracetam PK; one compartment, with
% absorption from virtual gut, and clearance from the central compartment.
% infusion parameter included for an intravenous delivery option
%

% note: y(1) is a concentration (mg/L), y(2) is an amount (mg)

ka = p(1);
kc = p(2);
Vd = p(3);
q = p(4);

dydt = zeros(3,1);    % initialize output vector

% 1 = Levetiracetam in central compartment (mg/L)
% 2 = Levetiracetam in gut virtual compartment (mg)
% 3 = Levetiracetam in cleared virtual compartment (mg)

dydt(1) = q/Vd + ka*y(2)/Vd - kc*y(1);
dydt(2) =      - ka*y(2);
dydt(3) =                   kc*y(1)*Vd;


