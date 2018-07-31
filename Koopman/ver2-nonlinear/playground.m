% just playing around

%% Showing that U' can flow the state forward in time
exnow = polyLift(xreal(1,:)');
fart(:,1) = exnow;
poop(:,1) = fart(2:3,1);
for i = 2:4
    fart(:,i) =  U' * exnow;
    xnow = fart(2:3,i);
    exnow = fart(:,i);
    
    poop(:,i) = xnow; 
end

figure
hold on
plot(poop(1,:))
plot(poop(2,:))
hold off