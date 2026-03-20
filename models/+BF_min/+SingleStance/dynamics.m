function f_S = dynamics(X, u, p)
    import BF_min.SingleStance.*;
    
    %% Initialization
    nq = length(X)/2;
    q = X(1:nq);       
    dq = X(nq+1:end);
    
    %% ddq 
    % cumpute ddq under minimal coordinate 
    M = massMatrix_m(X, p);  
    F = forces_m(X, u, p);   
    ddq = M \ F;
    
    %% dynamic
    f_S = [dq;ddq];
end


% function [f_S, lambda] = dynamics(X, u, p)
%     import BF_min.SigleStance.*;
%     %% initialization
%     % floating base
%     nq_f = length(X)/2;
%     q_f = X(1:nq_f);       
%     dq_f = X(nq_f+1:end);
% 
%     % mininal coordinate
%     q_m = q_f(3:end);       
%     dq_m = dq_f(3:end);
%     X_m = [q_m;dq_m];
% 
%     %% System matrix
%     M_m = massMatrix_m(X, p);  
%     F_m= force_m(X, u, p);
% 
%     M_f = massMatrix_f(X, p);  
%     F_f= force_f(X, u, p);   
% 
%     J_Stance = J_StanceContact(X, p);
%     dJ_Stance = dJ_StanceContact(X, p);
% 
%     % gamma = dJ* dq
%     gamma = dJ_Stance * dq;
% 
%     A = [M, -J_Stance';
%          J_Stance, zeros(size(J_Stance,1))];
%     B = [F_cg;
%          -gamma];
% 
%     sol = A \ B;
%     ddq = sol(1:nq);
%     lambda = sol(nq+1:end);
% 
%     %% dynamic
%     f_S = [dq;ddq];
% end

