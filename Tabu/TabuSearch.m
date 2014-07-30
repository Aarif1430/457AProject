%function  [BestSoln BestSolnCost] = TabuSearch( ...
%                ProbData, TabuLength, NumIterations, ...      
%                GenInitialSolnFn, GetBestNeighbourSolnFn, matrixSize, numOfTurbine)
function  [BestSoln BestSolnCost] = TabuSearch(TabuLength, NumIterations,matrixSize, numOfTurbine)


global size gridSize windVel rotorRadius N windSpeedMatrix
size= matrixSize;
gridSize = 80;
windVel=12;
rotorRadius=20;
N=numOfTurbine;
            
% This function implements the tabu search algorithm.
%
% Inputs:
%   ProbData: The data of the problem to be solved.
%   TabuLength: The length of the tabu list
%   NumIterations: The maximum number of iterations
%   GenInitialSolnFn: A handle to a function that generates an initial
%                     solution to the problem.
%   GetBestNeighbourSolnFn: A hanlde to a function that generates the 
%                         neighbourhood of a given solution and update
%                         the best neighborhood.
%
% Outputs:
%   BestSoln: The best solution obtained
%   BestSolnCost: The best solution cost

windSpeedMatrix = initWindSpeedMatrix(size);

% Generate the initial solution
[Soln SolnCost TabuList] = GenInitialSln(size,numOfTurbine);

% Set the best solution to the initial solution
BestSoln = Soln;
BestSolnCost = SolnCost;

for nIt = 1 : NumIterations
    % Get the best solution in the neighbourhood of the current solution
    % avoiding Tabu moves
    %[Soln SolnCost TabuList] = feval(GetBestNeighbourSolnFn, ...
    %                            Soln, TabuList, TabuLength, BestSolnCost);
    [Soln SolnCost TabuList] = GetBestNeighbourSolnFn(Soln, TabuList, TabuLength, BestSolnCost);
            
    % Update the best solution
    if SolnCost < BestSolnCost
        BestSoln = Soln;
        BestSolnCost = SolnCost;
    end
    %BestSoln
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

function [initialM, result, TabuList] = GenInitialSln(size, numTurbine)
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
    TabuList=zeros(size,size);
end

function [BestNeighbour,BestNeighbourCost,TabuList]=GetBestNeighbourSolnFn(Soln, TabuList, TabuLength, BestSolnCost) 
    bannedMoves=zeros(2, 1);
    mapsize=size(Soln,1);

    while true
        BestNeighbour = Soln;
        BestNeighbourCost = Inf;
        curBestCost = Inf;
        curBestSoln = Soln;
        potentialTabu = [0,0];
        
        for i = 1:mapsize
            for j = 1:mapsize
                if Soln(i,j)==1
                    %upper left
                    if(i-1>0 && j-1>0 && Soln(i-1,j-1)~=1 && checkBannedMoves(bannedMoves, i-1, j-1)==0)
                        matrix = Soln;

                        matrix(i-1,j-1)=1;
                        matrix(i,j)=0;

                        cost = CalculateCostFunc(matrix);
                        if(cost < curBestCost)
                            curBestCost = cost;
                            curBestSoln = matrix;
                            potentialTabu=[i-1,j-1];
                        end
                    end

                    %left
                    if(j-1>0 && Soln(i,j-1)~=1 && checkBannedMoves(bannedMoves, i, j-1)==0)
                        matrix = Soln;

                        matrix(i,j-1)=1;
                        matrix(i,j)=0;

                        cost = CalculateCostFunc(matrix);
                        if(cost < curBestCost)
                            curBestCost = cost;
                            curBestSoln = matrix;
                            potentialTabu=[i,j-1];
                        end
                    end

                    %lower left
                    if(i+1<=mapsize && j-1>0 && Soln(i+1,j-1)~=1 && checkBannedMoves(bannedMoves, i+1, j-1)==0)
                        matrix = Soln;

                        matrix(i+1,j-1)=1;
                        matrix(i,j)=0;

                        cost = CalculateCostFunc(matrix);
                        if(cost < curBestCost)
                            curBestCost = cost;
                            curBestSoln = matrix;
                            potentialTabu=[i+1,j-1];
                        end
                    end

                    %lower
                    if( i+1<=mapsize && Soln(i+1,j)~=1 && checkBannedMoves(bannedMoves, i+1, j)==0)
                        matrix = Soln;

                        matrix(i+1,j)=1;
                        matrix(i,j)=0;

                        cost = CalculateCostFunc(matrix);
                        if(cost < curBestCost)
                            curBestCost = cost;
                            curBestSoln = matrix;
                            potentialTabu=[i+1,j];
                        end
                    end

                    %lower right
                    if( i+1<=mapsize && j+1<=mapsize && Soln(i+1,j+1)~=1 && checkBannedMoves(bannedMoves, i+1, j+1)==0)
                        matrix = Soln;

                        matrix(i+1,j+1)=1;
                        matrix(i,j)=0;

                        cost = CalculateCostFunc(matrix);
                        if(cost < curBestCost)
                            curBestCost = cost;
                            curBestSoln = matrix;
                            potentialTabu=[i+1,j+1];
                        end
                    end

                    %right
                    if( j+1<=mapsize && Soln(i,j+1)~=1 && checkBannedMoves(bannedMoves, i, j+1)==0)
                        matrix = Soln;

                        matrix(i,j+1)=1;
                        matrix(i,j)=0;

                        cost = CalculateCostFunc(matrix);
                        if(cost < curBestCost)
                            curBestCost = cost;
                            curBestSoln = matrix;
                            potentialTabu=[i,j+1];
                        end
                    end

                    %upper right
                    if( i-1>0 && j+1<=mapsize && Soln(i-1,j+1)~=1 && checkBannedMoves(bannedMoves, i-1, j+1)==0)
                        matrix = Soln;

                        matrix(i-1,j+1)=1;
                        matrix(i,j)=0;

                        cost = CalculateCostFunc(matrix);
                        if(cost < curBestCost)
                            curBestCost = cost;
                            curBestSoln = matrix;
                            potentialTabu=[i-1,j+1];
                        end
                    end

                    %upper
                    if( i-1>0 && Soln(i-1,j)~=1 && checkBannedMoves(bannedMoves, i-1, j)==0)
                        matrix = Soln;

                        matrix(i-1,j)=1;
                        matrix(i,j)=0;

                        cost = CalculateCostFunc(matrix);
                        if(cost < curBestCost)
                            curBestCost = cost;
                            curBestSoln = matrix;
                            potentialTabu=[i-1,j];
                        end
                    end

                end%end for if Soln(i,j)==1
            end%end for j loop
        end%end for i loop

        %update other tabu value
        TabuList=TabuList-1;

        %aspiration criterion
        if or(TabuList(potentialTabu(1),potentialTabu(2))<=0, curBestSoln < BestSolnCost)
            BestNeighbour = curBestSoln;
            BestNeighbourCost = curBestCost;
            if potentialTabu(1)~=0
                TabuList(potentialTabu(1),potentialTabu(2))=TabuLength;
            end
            break;
        else
            index=size(bannedMoves,2); %num of column
            index=index+1;
            bannedMoves(1, index) = potentialTabu(1);
            bannedMoves(2, index) = potentialTabu(2);
        end
    end
end

function result = checkBannedMoves(bannedMoves, i, j)
    result = 0;

    for iter=1:size(bannedMoves,2)
        if i==bannedMoves(1,iter) && j==bannedMoves(2,iter)
            result=1;
            break;
        end
    end
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