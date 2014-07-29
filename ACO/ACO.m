% Ant Colony
%
% 

function [curBestSolnCost, curBestSol] = ACO(iterations, matrixSize, numOfTurbine, numOfAnt)  
  global size gridSize windVel rotorRadius N pheromoneMatrix alpha beta r0 windSpeedMatrix
  size= matrixSize;
  gridSize = 80;
  windVel=12;
  rotorRadius=20;
  N=numOfTurbine; 
  alpha=1;
  beta=1;
  r0=0.5;
  rho1=0.3;
  rho2=1;

  windSpeedMatrix = initWindSpeedMatrix(size);
  pheromoneMatrix=ones(size);%initial pheromone concentration is 1
  %pheromoneMatrix=pheromoneMatrix*100;
  curBestSol=zeros(size);
  curBestSolnCost=Inf;

  for it = 1:iterations
    for k = 1:numOfAnt
        [newSol, newSolCost] = GenSln(size, N);
        if newSolCost < curBestSolnCost
            curBestSol = newSol;
            curBestSolnCost = newSolCost;
        end
    end % end for ants

    %update pheromone concentrate
    reinforce=newSol.*pheromoneMatrix;
    reinforce=(1-rho1)*reinforce+newSol*rho2*(0.000001/curBestSolnCost);

    decay=(1-newSol).*pheromoneMatrix;
    decay=(1-rho1)*decay;

    pheromoneMatrix=reinforce+decay;

  end % end iterations
  
  
  
end

% initialize the wind park with different wind speeds
function windSpeedMatrix = initWindSpeedMatrix(size)
    global windVel
    % init a N by N matrix
    m = ones(size);
    m = m*12;
    
    for i=1:floor(size/4)
        for j=0:3
            windDiff = 2 * j;
            m(i+j,:) = windDiff + windVel;
        end
    end
    windSpeedMatrix=m;
end

% generate solution based on pheromone concentration
function [m, cost] = GenSln(size, numTurbine)
    global pheromoneMatrix alpha beta r0 windSpeedMatrix

    probVector=zeros(1,size*size);
    m=zeros(size,size);
    count=0;
    resetIndex=0;

    %calculate probability function for each location in the wind park
    for i = 1:size
        for j = 1:size
            desirability = windSpeedMatrix(i, j);
            probVector(1, (i-1)*size+j) = (pheromoneMatrix(i, j)^alpha) * ((desirability)^beta);
        end
    end
    [ma i]=min(probVector);
    %ma
    
    %probVector

    while count<numTurbine
        if rand() >= r0
            %Roulette Wheel method
            probVector = probVector / sum(probVector);
            cumulativeProbVector = cumsum(probVector);
            randnum = rand();
            for ii = 0:length(cumulativeProbVector)
                if ii == 0
                    if 0 < randnum && cumulativeProbVector(ii+1) >= randnum
                        loc = 1;
                        resetIndex = 1;
                    end
                elseif (cumulativeProbVector(ii) < randnum) && (cumulativeProbVector(ii+1) >= randnum)
                    loc = ii+1;
                    resetIndex = ii+1;
                end 
            end
        else
            [value loc] = max(probVector);
            resetIndex = loc;
        end

        %place a turbine at this chosen location in the wind park
        %first check if the location is already taken
        i = ceil(loc/size);
        j = mod(loc,size);

        if j==0
            j = size;
            if (i~=1)
                i = i - 1;
            end
        end
        if m(i,j)==0 
            m(i,j)=1;
            count=count+1;
            % set the probability of the recently chosen location to 0 to avoid choosing it again
            probVector(resetIndex) = 0;
        end
    end

    cost = CalculateCostFunc(m);

end


function result = CalculateCostFunc(m)
    N=sum(sum(m));
    cost = N*(2/3+1/3*exp(-0.00174*N^2));
    %cost
    result = cost / CalculateTotalPower(m);
end

function pwr = CalculateTotalPower(m)
    global size
    
    totalpower=0;
    for i = 1:size
        for j = 1:size
            if(m(i,j)==1)
                totalpower = totalpower+CalculateSingleTurbinePower(calculate_velocity(m,i,j));
            end
        end
    end
    pwr=totalpower*31536000;
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

function vel = calculate_velocity(matrix, x, y) 
    global gridSize size windVel windSpeedMatrix

    %thrus coefficient of the turbine
    ct = 0.88;
    k = 2;
    R = 20;
    
    vel_def_total = 0;
    vel = windSpeedMatrix(x, y);
    
    for i = 1 : size
        for j = 1 : y-1
            if matrix(i,j)==1
                if check_wake(i*gridSize, j*gridSize, x*gridSize, y*gridSize)==1
                    vel_def_cur = (1-sqrt(1-ct))/(1+k*(y-j)*gridSize/R)^2;
                    vel_def_total = vel_def_total+vel_def_cur^2;
                end
            end
        end
    end
    vel_def = sqrt(vel_def_total);
    % do not update velocity if turbine is not affected by wake loss
    if (vel_def ~= 0)
        vel = vel * (1-vel_def);
    end
end



%check if turbine i(x2,y2) is in wake of turbine j(x1,y1)
function check = check_wake(y1,x1,y2,x2)
    %wind direction is 0 degree
    theta =0;
    alpha = atan(2);
    R=20;
    k=2;
    
    numerator = (x2-x1)*cos(theta)+(y2-y1)*sin(theta)+R/k;
    denominator = sqrt( (x2-x1+R/k*cos(theta))^2 + (y2-y1+R/k*cos(theta))^2);
    beta = acos( numerator/denominator );

    if(beta<alpha)
        check=1;
    else
        check=0;
    end
end
