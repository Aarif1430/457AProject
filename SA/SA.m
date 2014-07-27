% Simulated Annealing Algorithm
%
% costs: history of cost for each iteration's schedule.
% bestSol: best solution for the given sT and fT.
%
% iterations: number of iteration for each T.
% sT: start T.
% fT: final T.
% alpha: alpha in geometric cooling schedule.



function [bestSolCost, bestSol] = SA(iterations, sT, fT, alpha, matrixSize, numOfTurbine)
                           
  global size gridSize windVel rotorRadius N windSpeedMatrix
  size= matrixSize;
  gridSize = 80;
  windVel=12;
  rotorRadius=20;
  N=numOfTurbine; 
  
  costsEnd = 0;
  costs = [];
  
  windSpeedMatrix = initWindSpeedMatrix(size);
  
  %generate a random initial soln
  [bestSol,bestSolCost] = GenInitialSln(size, N);
  acceptedSol = bestSol;
  acceptedSolCost = bestSolCost;
  T = sT;
  while T > fT
    for i=1:iterations
      [newSol, newSolCost] = GetBestNeighbourSolnFn(acceptedSol);
      deltaCost = newSolCost - acceptedSolCost;
      if deltaCost < 0
        acceptedSol = newSol;
        acceptedSolCost = newSolCost;
      else
        randVal = rand(1);
        p = exp(-1*deltaCost / T);
        if p > randVal
          acceptedSol = newSol;
          acceptedSolCost = newSolCost;
        end
      end

      % record the cost value in to history
      costsEnd = costsEnd + 1;
      costs(costsEnd) = acceptedSolCost;

      % Update current best value
      if acceptedSolCost < bestSolCost
         bestSol = acceptedSol;
         bestSolCost = acceptedSolCost;
      end
    end
    T = T * alpha; % cooling
  end
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

function [BestNeighbour,BestNeighbourCost]=GetBestNeighbourSolnFn(Soln)
    global size
    
    BestNeighbour = Soln;
    BestNeighbourCost = Inf;
    curBestCost = Inf;
    curBestSoln = Soln;
    for i = 1:size
        for j = 1:size
            if Soln(i,j)==1
                %upper left
                if(i-1>0 && j-1>0 && Soln(i-1,j-1)~=1)
                    matrix = Soln;
                    
                    matrix(i-1,j-1)=1;
                    matrix(i,j)=0;
                    
                    cost = CalculateCostFunc(matrix);
                    if(cost < curBestCost)
                        curBestCost = cost;
                        curBestSoln = matrix;
                    end
                end
                
                %left
                if(j-1>0 && Soln(i,j-1)~=1)
                    matrix = Soln;
                    
                    matrix(i,j-1)=1;
                    matrix(i,j)=0;
                    
                    cost = CalculateCostFunc(matrix);
                    if(cost < curBestCost)
                        curBestCost = cost;
                        curBestSoln = matrix;
                    end
                end
                
                %lower left
                if(i+1<=size && j-1>0 && Soln(i+1,j-1)~=1)
                    matrix = Soln;
                    
                    matrix(i+1,j-1)=1;
                    matrix(i,j)=0;
                    
                    cost = CalculateCostFunc(matrix);
                    if(cost < curBestCost)
                        curBestCost = cost;
                        curBestSoln = matrix;
                    end
                end
                
                %lower
                if( i+1<=size && Soln(i+1,j)~=1)
                    matrix = Soln;
                    
                    matrix(i+1,j)=1;
                    matrix(i,j)=0;
                    
                    cost = CalculateCostFunc(matrix);
                    if(cost < curBestCost)
                        curBestCost = cost;
                        curBestSoln = matrix;
                    end
                end
                
                %lower right
                if( i+1<=size && j+1<=size && Soln(i+1,j+1)~=1)
                    matrix = Soln;
                    
                    matrix(i+1,j+1)=1;
                    matrix(i,j)=0;
                    
                    cost = CalculateCostFunc(matrix);
                    if(cost < curBestCost)
                        curBestCost = cost;
                        curBestSoln = matrix;
                    end
                end
                
                %right
                if( j+1<=size && Soln(i,j+1)~=1)
                    matrix = Soln;
                    
                    matrix(i,j+1)=1;
                    matrix(i,j)=0;
                    
                    cost = CalculateCostFunc(matrix);
                    if(cost < curBestCost)
                        curBestCost = cost;
                        curBestSoln = matrix;
                    end
                end
                
                %upper right
                if( i-1>0 && j+1<=size && Soln(i-1,j+1)~=1)
                    matrix = Soln;
                    
                    matrix(i-1,j+1)=1;
                    matrix(i,j)=0;
                    
                    cost = CalculateCostFunc(matrix);
                    if(cost < curBestCost)
                        curBestCost = cost;
                        curBestSoln = matrix;
                    end
                end
                
                %upper
                if( i-1>0 && Soln(i-1,j)~=1)
                    matrix = Soln;
                    
                    matrix(i-1,j)=1;
                    matrix(i,j)=0;
                    
                    cost = CalculateCostFunc(matrix);
                    if(cost < curBestCost)
                        curBestCost = cost;
                        curBestSoln = matrix;
                    end
                end
          
            end%end for if Soln(i,j)==1
        end%end for j loop
    end%end for i loop

    BestNeighbour = curBestSoln;
    BestNeighbourCost = curBestCost;
    
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

function vel = calculate_velocity(matrix, x, y) 
    global gridSize size windVel windSpeedMatrix

    %thrus coefficient of the turbine
    ct = 0.88;
    k = 2;
    R = 20;
    
    vel_def_total = 0;
    vel = windSpeedMatrix(x,y);
    
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
