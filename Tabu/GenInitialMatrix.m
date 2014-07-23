function GetInitialMatrix()
    
end

function [initialM, result] = GenInitialSln(size, numTurbine)
    m=zeros(size,size);
    count=0;
    while count<numTurbine
        i=ceil(rand()*size);
        j=ceil(rand()*size);
        if( m(i,j)==0 )
            m(i,j)=1;
            count=count+1;
        end
    end
    
    initialM = m;
    result = CalculateCostFunc(m);
end

function result = CalculateCostFunc(m)
    N=sum(sum(m));
    cost = N*(2/3+1/3*exp(-0.00174*N^2));
    result = cost / CalculateTotalPower(m);
end

function pwr = CalculateTotalPower(m)
    global matrixSize
    
    totalpower=0;
    for i = 1:matrixSize
        for j = 1:matrixSize
            if(m(i,j)==1)
                totalpower = totalpower+CalculateSingleTurbinePower(calculate_velocity(m,i,j));
            end
        end
    end
    pwr=totalpower;
end

function pwr = CalculateSingleTurbinePower(vel)
    global rotorRadius
    A = pi*rotorRadius^2;
    rho = 1.2;
    Cp = 0.35;
    Ng = 0.7;
    Nb = 0.95;
    
    pwr = 0.5 * rho * A * Cp * vel^3 * Ng * Nb;
end
