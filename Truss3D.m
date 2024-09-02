clear all
clc
format long

%%% FEA OF 3D TRUSS

%%% Data
n = 48; % Number of truss elements
modElasticity = 30*10^6*ones(n,1); % modulus of elasticity
 area =  6*ones(n,1); %[25 12 25 12 1 4 17 17 17 5]'; % Area of truss elements column vector
coOrd = [0 0 0;
         117.6 0 0;
         0 0 -117.6;
         117.6 0 -117.6;
         0 132 0;
         117.6 132 0;
         0 132 -117.6;
         117.6 132 -117.6;
         0 288 0;
         117.6 288 0;
         0 288 -117.6;
         117.6 288 -117.6;
         0 444 0;
         117.6 444 0;
         0 444 -117.6;
         117.6 444 -117.6];% Coordinates of nodes (xOrd yOrd zord) in node order 1,2,3 etc
       
nodeNo = [1:16]; % Node number 

for i=1:length(nodeNo) % DOF associated with the node number
DOF(i,:) = [3*nodeNo(i)-2 3*nodeNo(i)-1 3*nodeNo(i)];
end

eleConnect = [5 6;
              6 8;
              8 7;
              7 5;
              9 10;
              10 12;
              12 11;
              11 9;
              13 14;
              14 16;
              16 15;
              15 13;
              1 6;2 5;2 8;4 6;4 7;3 8;1 7;3 5;5 10;6 9;7 12;11 8;5 11;9 7;12 6;10 8;9 14;10 13;15 12;11 16;10 16;12 14;9 15;11 13;
              1 5;2 6;3 7;4 8;5 9;6 10;7 11;8 12;9 13;10 14;11 15;12 16];
              
              % Element connectivity element rowwise 1,2,3 etc
              % First column(local node 1)
              % Second column (local node 2)
             

for i = 1:n % Length of element sqrt((X2-X1)^2+(Y2-Y1)^2)
    le(i,1) = sqrt((coOrd(eleConnect(i,2),1)-coOrd(eleConnect(i,1),1))^2+(coOrd(eleConnect(i,2),2)-coOrd(eleConnect(i,1),2))^2+(coOrd(eleConnect(i,2),3)-coOrd(eleConnect(i,1),3))^2);
end

for i = 1:n % Direction cosine l = (X2-X1)/le  m = (Y2-Y1)/le  n = (Z2-Z1)/le  element rowwise
    dirCosines(i,:) = [(coOrd(eleConnect(i,2),1)-coOrd(eleConnect(i,1),1))/le(i) (coOrd(eleConnect(i,2),2)-coOrd(eleConnect(i,1),2))/le(i) (coOrd(eleConnect(i,2),3)-coOrd(eleConnect(i,1),3))/le(i)];
end

% Stiffness matrix
for k = 1:n % element stiffness matrix
eleStiff(:,:,k) = (modElasticity(k)*area(k)/le(k)) * [dirCosines(k,1)^2 dirCosines(k,1)*dirCosines(k,2) dirCosines(k,1)*dirCosines(k,3)  -dirCosines(k,1)^2 -dirCosines(k,1)*dirCosines(k,2) -dirCosines(k,1)*dirCosines(k,3);
                                                      dirCosines(k,1)*dirCosines(k,2) dirCosines(k,2)^2 dirCosines(k,2)*dirCosines(k,3) -dirCosines(k,1)*dirCosines(k,2) -dirCosines(k,2)^2 -dirCosines(k,2)*dirCosines(k,3);
                                                      dirCosines(k,1)*dirCosines(k,3) dirCosines(k,2)*dirCosines(k,3) dirCosines(k,3)^2 -dirCosines(k,1)*dirCosines(k,3) -dirCosines(k,2)*dirCosines(k,3) -dirCosines(k,3)^2;
                                                      -dirCosines(k,1)^2 -dirCosines(k,1)*dirCosines(k,2) -dirCosines(k,1)*dirCosines(k,3) dirCosines(k,1)^2 dirCosines(k,1)*dirCosines(k,2) dirCosines(k,1)*dirCosines(k,3);
                                                     -dirCosines(k,1)*dirCosines(k,2) -dirCosines(k,2)^2 -dirCosines(k,2)*dirCosines(k,3)  dirCosines(k,1)*dirCosines(k,2) dirCosines(k,2)^2 -dirCosines(k,2)*dirCosines(k,3);
                                                     -dirCosines(k,1)*dirCosines(k,3) -dirCosines(k,2)*dirCosines(k,3) -dirCosines(k,3)^2 dirCosines(k,1)*dirCosines(k,3) dirCosines(k,2)*dirCosines(k,3) dirCosines(k,3)^2];
end

for k = 1:length(le)
           nodeLocal1 = eleConnect(k,1);
           nodeLocal2 = eleConnect(k,2);
           dof1 = DOF(nodeLocal1,:);  % dof associated with nodes
           dof2 = DOF(nodeLocal2,:);
           assStiff(:,:,k) = zeros(max(max(DOF))); % Assemble Global Stiffness Matrix
           assStiff(dof1(1),dof1(1),k) = eleStiff(1,1,k);
           assStiff(dof1(1),dof1(2),k) = eleStiff(1,2,k);
           assStiff(dof1(1),dof1(3),k) = eleStiff(1,3,k);
           assStiff(dof1(1),dof2(1),k) = eleStiff(1,4,k);
           assStiff(dof1(1),dof2(2),k) = eleStiff(1,5,k); 
           assStiff(dof1(1),dof2(3),k) = eleStiff(1,6,k); 

           assStiff(dof1(2),dof1(1),k) = eleStiff(2,1,k);
           assStiff(dof1(2),dof1(2),k) = eleStiff(2,2,k);
           assStiff(dof1(2),dof1(3),k) = eleStiff(2,3,k);
           assStiff(dof1(2),dof2(1),k) = eleStiff(2,4,k);
           assStiff(dof1(2),dof2(2),k) = eleStiff(2,5,k); 
           assStiff(dof1(2),dof2(3),k) = eleStiff(2,6,k); 

           assStiff(dof1(3),dof1(1),k) = eleStiff(3,1,k);
           assStiff(dof1(3),dof1(2),k) = eleStiff(3,2,k);
           assStiff(dof1(3),dof1(3),k) = eleStiff(3,3,k);
           assStiff(dof1(3),dof2(1),k) = eleStiff(3,4,k);
           assStiff(dof1(3),dof2(2),k) = eleStiff(3,5,k); 
           assStiff(dof1(3),dof2(3),k) = eleStiff(3,6,k); 

           assStiff(dof2(1),dof1(1),k) = eleStiff(4,1,k);
           assStiff(dof2(1),dof1(2),k) = eleStiff(4,2,k);
           assStiff(dof2(1),dof1(3),k) = eleStiff(4,3,k);
           assStiff(dof2(1),dof2(1),k) = eleStiff(4,4,k);
           assStiff(dof2(1),dof2(2),k) = eleStiff(4,5,k); 
           assStiff(dof2(1),dof2(3),k) = eleStiff(4,6,k);
           
           assStiff(dof2(2),dof1(1),k) = eleStiff(5,1,k);
           assStiff(dof2(2),dof1(2),k) = eleStiff(5,2,k);
           assStiff(dof2(2),dof1(3),k) = eleStiff(5,3,k);
           assStiff(dof2(2),dof2(1),k) = eleStiff(5,4,k);
           assStiff(dof2(2),dof2(2),k) = eleStiff(5,5,k); 
           assStiff(dof2(2),dof2(3),k) = eleStiff(5,6,k);

           assStiff(dof2(3),dof1(1),k) = eleStiff(6,1,k);
           assStiff(dof2(3),dof1(2),k) = eleStiff(6,2,k);
           assStiff(dof2(3),dof1(3),k) = eleStiff(6,3,k);
           assStiff(dof2(3),dof2(1),k) = eleStiff(6,4,k);
           assStiff(dof2(3),dof2(2),k) = eleStiff(6,5,k); 
           assStiff(dof2(3),dof2(3),k) = eleStiff(6,6,k);
end

globStiff=  assStiff(:,:,1);
for i = 2:length(le)
    globStiff = globStiff + assStiff(:,:,i); % Global stiiffness matrix
end

finalStiffMat = globStiff;

disp('Global stiffness matrix is');
disp(globStiff);

%%% Joint load Vector
P = zeros(max(max(DOF)),1);% Joint load acting on the truss
jointLoadDof = [37 43];
jointLoad = [5280 5280];
for i = 1:length(jointLoadDof)
    P(jointLoadDof(i)) = jointLoad(i);
end

%%% Temperature/lack of fit load vector
switchTemp = 0; % If temperature or lack of fit is considered put switch=1 
    if switchTemp == 1
        membaffected = [1]; % Member effected by temperature or lack of fit
        alpha = .01; % Coefficient of thermal expansion
        deltaT = 1 ; % Change in temperature
        initialstrain = alpha*deltaT*[1 0]; % if lack of fit put value and if temoperature effect put alpha*deltaT member row wise 1,2,3 etc
    for k = 1:length(membaffected) % element temperature/lack of fit vector
        elemTempLoad(:,:,k) = modElasticity(membaffected(k))*area(membaffected(k))*initialstrain(membaffected(k))*[-dirCosines(membaffected(k),1);
                                                                                                  -dirCosines(membaffected(k),2);
                                                                                                  -dirCosines(membaffected(k),3);
                                                                                                   dirCosines(membaffected(k),1);
                                                                                                   dirCosines(membaffected(k),2);
                                                                                                   dirCosines(membaffected(k),3)
                                                                                                   ];
    end

for k = 1:length(membaffected)
            assForce(:,:,k) = zeros(max(max(DOF)),1); % Assemble Global Force Vector
            nodeLocal1 = eleConnect(membaffected(k),1);
            nodeLocal2 = eleConnect(membaffected(k),2);
            dof1 = DOF(nodeLocal1,:);  % dof associated with nodes
            dof2 = DOF(nodeLocal2,:);
            assForce(dof1(1),1,k) = elemTempLoad(1,1,k);
            assForce(dof1(2),1,k) = elemTempLoad(2,1,k);
            assForce(dof1(3),1,k) = elemTempLoad(3,1,k);
            assForce(dof2(1),1,k) = elemTempLoad(4,1,k);
            assForce(dof2(2),1,k) = elemTempLoad(5,1,k);
            assForce(dof2(3),1,k) = elemTempLoad(6,1,k);

end
     
globForce=  assForce(:,:,1);
for i = 2:length(membaffected)
    globForce = globForce + assForce(:,:,i);
end
    else
globForce = zeros(max(max(DOF)),1);
initialstrain = zeros(length(le),1);
    end
globForce = globForce + P; % GLobal force vector
finalForceVec = globForce;

disp('Global Force vector is');
disp(globForce)

format longg

%%% Constriants
% Constraint type b1Qp1 + b2Qp2 = bo
b1b2 = [1;1;1;1;1;1;1;1;1;1;1;1]; % coefficients of contraint (number of rows indicates a constraint equations)
bo = zeros(12,1);
constDOF = [1;2;3;4;5;6;7;8;9;10;11;12
    ]; % DOF constrained (one row one equation)
[row,cols] = size(constDOF);
C = max(max(globStiff))*10^4; % Penalty method

% modification to stiffness and force matrix
for i = 1:row
    for j = 1:cols;
        globStiff(constDOF(i,j),constDOF(i,j)) = globStiff(constDOF(i,j),constDOF(i,j))+C*b1b2(i,j)*b1b2(i,j);
    end
        if cols > 1
        globStiff(constDOF(i,j-1),constDOF(i,j)) = globStiff(constDOF(i,j-1),constDOF(i,j))+C*b1b2(i,j-1)*b1b2(i,j);
        globStiff(constDOF(i,j),constDOF(i,j-1)) = globStiff(constDOF(i,j-1),constDOF(i,j));
        end
end
      
for i = 1:row
    for j = 1:cols
    globForce(constDOF(i,j)) = globForce(constDOF(i,j))+C*bo(i)*b1b2(i,j);
    end
end


% solve equation

n = length(globForce); % Number of equation

% Forward elimination (reduction of A,b)
for k = 1:n-1
    for i = k+1:n
        c = globStiff(i,k)/globStiff(k,k);
        for j=k:n
           globStiff(i,j) = globStiff(i,j)-c*globStiff(k,j);
        end
        globForce(i) = globForce(i)-c*globForce(k);
    end
end


% Back substitution
globForce(n) = globForce(n)/globStiff(n,n);
for ii = 1:n-1
    i = n-ii;
    sum = 0;
    for j = i+1:n ;
        sum = sum + globStiff(i,j)*globForce(j);
    end
    globForce(i) = (globForce(i)-sum)/globStiff(i,i);
end
finalDispVec = globForce;
disp('The displacements are')
disp(globForce); % Vector giving displacement


% Stress calculation
for i = 1:length(le)
stress(i,1) = modElasticity(i)*((1/le(i))*[-dirCosines(i,1) -dirCosines(i,2) -dirCosines(i,3) dirCosines(i,1) dirCosines(i,2) dirCosines(i,3)]*[globForce(DOF(eleConnect(i,1),1)); globForce(DOF(eleConnect(i,1),2)); globForce(DOF(eleConnect(i,1),3)); globForce(DOF(eleConnect(i,2),1)); globForce(DOF(eleConnect(i,2),2)); globForce(DOF(eleConnect(i,2),3))]-initialstrain(i)); 
end

disp('Stress in elements are')
disp(stress)

% Force calculation
forceMember = stress.*area; % Force in member tension (+) compression (-)
disp('Force in members')
disp(forceMember)

% Reaction calculation R = KQ-F (only constrained DOF)
for i = 1:row
    for j = 1:cols
        k = constDOF(i,j);
reaction(i,1) = finalStiffMat(k,:)*finalDispVec - finalForceVec(k);
    end
end
disp('Reactions are')
disp(reaction)
           




