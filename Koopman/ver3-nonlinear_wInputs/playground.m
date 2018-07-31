% just playing around

%%
exnow = polyLift(xreal(1,:)');
for i = 1:length(treal)
    fart(:,i) = blkdiag(0, eye(2), zeros(17)) * U' * exnow;
    xnow = fart(:,i);
    exnow = polyLift(xnow);
end

figure
plot(treal, fart)
    