function F = P_dot(A,B,Q,R,P)

F = -(-A'*P-P*A+P*B*inv(R)*B'*P-Q);

end