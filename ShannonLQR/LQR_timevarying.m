function K = LQR_timevarying(A,B,Q,R,t,dt)
% This function computes a time-varying gain matrix K(t)
% given cells A{}, B{}, Q{}, and R{} and a time vector
%store each element of the gain matrix in a cell K{1}, K{2}, etc

% 10/19/2016
%
% first, solve for P(t) using Heun's method
P=cell(length(t),1); P_intermed=cell(length(t),1); K=cell(length(t),1);
P{1}=Q{1};
%
for i=1:length(t)-1
    %initial prediction P_intermed made with Euler's method (not very
    %accurate)
    P_intermed{i+1}=P{i}+dt*P_dot(A{i},B{i},Q{i},R{i},P{i});
    %final approximation of P{i+1} is achieved by using an average of the
    %two slopes from P_i and P_intermed
    P{i+1}=P{i}+dt/2*(P_dot(A{i},B{i},Q{i},R{i},P{i})+...
        P_dot(A{i+1},B{i+1},Q{i+1},R{i+1},P_intermed{i+1}));
    %Now determine K{i}
    K{i}=inv(R{i})*B{i}'*P{i};
end
%K{end} isn't defined in the loop (only goes to length(t)-1)
K{end}=inv(R{end})*B{end}'*P{end};
end