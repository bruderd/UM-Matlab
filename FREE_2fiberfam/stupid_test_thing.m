% Filter out cases where gama is not bigger in magnitude
qplot = zeros(181,3);
count = 1;
for k = 1:length(qual(:,1))
    
    count = count + 1;
    if (abs(qual(k,1)) > abs(qual(k,2)))
        qplot(count,:) = [qual(k,1), qual(k,2), 200*qual(k,3) + 100*qual(k,4)];
    end

%     qplot(count,:) = [qual(k,1), qual(k,2), 200*qual(k,3) + 100*qual(k,4)];

end