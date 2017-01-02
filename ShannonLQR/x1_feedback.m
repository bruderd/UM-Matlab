function f = x1_feedback(A,B,K,x1,x2)

G=A-B*K;
f = G(1,1)*x1+G(1,2)*x2;

end