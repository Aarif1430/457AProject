x=[];
y=[];
count=1;
for numTurbine=5:20
    %run each case 5 times and take average
    sum=0;
    for i=1:5
        [bestSoln, solnCost] = TabuSearch(5, 100 ,100,numTurbine );
        sum=solnCost+sum;
    end
    x(count) = numTurbine;
    y(count) = sum/5; %take average
    count=count+1;
end

plot(x, y);