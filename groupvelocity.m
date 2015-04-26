pas=[435 335 315 320 350];
nee=[550 430 405 407 460];
U=331./(nee-pas);
T=[20 30 40 50 60];
figure;
plot(T(1:4),U(1:4),'r*-',T(1:4),[3.62 3.88 3.72 3.92],'b*-');
title('U,c--T');
xlabel('period T');
ylabel('velocity');
legend('group velocity','phase velocity')