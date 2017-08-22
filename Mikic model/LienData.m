%% Lien's experiment data
% B2, B5, B7
% P         1262,   12589.1, 38660.3
% Delta_T   15.74,  10.67,  9.01
% Ja        2690.0, 201.1,  57.5
% Tinf      299.3,  334.24,  357.23

%B2
t2 = [0.0004 0.0009 0.0014 0.0019 0.0024 ...
      0.0029 0.0034 0.0039 0.0044 0.0049 ...
      0.0054 0.0059 0.0064 0.0069 0.0074 ...
      0.0079];
R2 = [0.0430 0.0974 0.1481 0.1963 0.2481 ...
      0.2979 0.3469 0.3963 0.4464 0.4968 ...
      0.5482 0.5993 0.6559 0.7153 0.7686 ...
      0.8174]*1e-2;
  
%B5
t5 = [0.0003 0.0008 0.0013 0.0018 0.0023 ...
      0.0028 0.0033 0.0038 0.0043 0.0048 ...
      0.0053 0.0058 0.0063 0.0068 0.0073 ...
      0.0078];
R5 = [0.0578 0.1446 0.2190 0.2904 0.3542 ...
      0.4155 0.4761 0.5317 0.5868 0.6370 ...
      0.6872 0.7326 0.7803 0.8366 0.8819 ...
      0.9257]*1e-2;
  
%B7
t7 = t2 + 0.00008;
R7 = [0.0626 0.1192 0.1572 0.1880 0.2142 ...
      0.2371 0.2569 0.2760 0.2943 0.3103 ...
      0.3271 0.3406 0.3561 0.3699 0.3862 ...
      0.3987]*1e-2;
  
figure;plot(t,Examination(:,1)*100,t,Examination(:,2)*100,t,Examination(:,3)*100)
hold on
scatter(t2, R2*100)
scatter(t5, R5*100)
scatter(t7, R7*100)
legend('B2 - Model','B5 - Model','B7 - Model','B2 - Experiment','B5 - Experiment','B7 - Experiment')
xlabel('time(s)');ylabel('R(cm)');
title('Radius development of gas bubble');