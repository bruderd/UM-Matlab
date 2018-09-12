function some_snapshotPairs = get_randsnapshotPairs( num , snapshotPairs )
%get_randsnapshotPairs: Randomly extracts num snapshot pairs
%   Detailed explanation goes here

some_snapshotPairs = struct;

totalPairs = length(snapshotPairs.x);

index = randi(totalPairs, num, 1);

some_snapshotPairs.x = snapshotPairs.x(index,:);
some_snapshotPairs.y = snapshotPairs.y(index,:);

end

