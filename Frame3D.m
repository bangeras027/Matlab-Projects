clear all
clc
format shortEng

%%% FEA OF 3D FRAME

%%% Data
n = 8; % Number of frame elements
modElasticity = 30*10^6*ones(n,1); % modulus of elasticity
poisRatio = .3; % possion ratio
shearModRig = modElasticity./(2*(1+poisRatio)); % modulus of rigidity
area = [6;6;6;6;3;3;3;3]; % Area of frame elements column vector
momInertiaZ = [51;51;51;51;17;17;17;17]; % moment of inertia about local z axis
momInertiaY = [3.75;3.75;3.75;3.75;1.26;1.26;1.26;1.26]; % moment of inertia about local y axis
polMomInertia = [.24;.24;.24;.24;.08;.08;.08;.08]; % Polar moment of inertia about local x axis

coOrd = [0 0 180;    % Coordinates of nodes (xOrd yOrd) in node order 1,2,3 etc
         0 180 180;
         120 180 180;
         120 0 180;
         0 0 0;
         0 180 0;
         120 180 0;
         120 0 0];

referenceNode = [60 90 90 ]; % Reference node should not lie in line connecting any node

nodeNo = [1:8]; % Node number 

for i=1:length(nodeNo) % DOF associated with the node number
DOF(i,:) = [6*nodeNo(i)-5 6*nodeNo(i)-4 6*nodeNo(i)-3 6*nodeNo(i)-2 6*nodeNo(i)-1 6*nodeNo(i)]; % [Translation Rotation]
end

eleConnect = [1 2; % Element connectivity element rowwise 1,2,3 etc First column(local node 1)  Second column (local node 2)
              4 3;
              5 6;
              8 7;
              2 3;
              2 6;
              6 7;
              3 7
             ];

for i = 1:n
DOFConnect(i,:) = [DOF(eleConnect(i,1),:) DOF(eleConnect(i,2),:)]; % DOF connectivity element wise
end

for i = 1:n % Length of element sqrt((X2-X1)^2+(Y2-Y1)^2)
    le(i,1) = sqrt((coOrd(eleConnect(i,2),1)-coOrd(eleConnect(i,1),1))^2+(coOrd(eleConnect(i,2),2)-coOrd(eleConnect(i,1),2))^2+(coOrd(eleConnect(i,2),3)-coOrd(eleConnect(i,1),3))^2);
end

% Calculation of direction cosines
for i = 1:n % Direction cosine l1 = (X2-X1)/le  m2 = (Y2-Y1)/le n2 = (Z2-Z1)/le  element rowwise
    dirCosinesx1(i,:) = [(coOrd(eleConnect(i,2),1)-coOrd(eleConnect(i,1),1))/le(i) (coOrd(eleConnect(i,2),2)-coOrd(eleConnect(i,1),2))/le(i) (coOrd(eleConnect(i,2),3)-coOrd(eleConnect(i,1),3))/le(i)];
end

for i = 1:n 
    l13(i,1) = sqrt((referenceNode(1)-coOrd(eleConnect(i,1),1))^2+(referenceNode(2)-coOrd(eleConnect(i,1),2))^2+(referenceNode(3)-coOrd(eleConnect(i,1),3))^2);
end

for i = 1:n % Direction cosine l1 = (X2-X1)/le  m2 = (Y2-Y1)/le n2 = (Z2-Z1)/le  element rowwise
    dirCosines13(i,:) = [(referenceNode(1) -coOrd(eleConnect(i,1),1))/l13(i) (referenceNode(2)-coOrd(eleConnect(i,1),2))/l13(i) (referenceNode(3)-coOrd(eleConnect(i,1),3))/l13(i)];
end
 
for i = 1:n 
dirCosinesx1XdirCosines13(i,:) = cross(dirCosinesx1(i,:),dirCosines13(i,:)); 
dirCosinesz1(i,:) = dirCosinesx1XdirCosines13(i,:) ./ sqrt(dirCosinesx1XdirCosines13(i,1)^2 + dirCosinesx1XdirCosines13(i,2)^2 + dirCosinesx1XdirCosines13(i,3)^2);
end

for i = 1:n
    dirCosinesy1(i,:) = cross(dirCosinesz1(i,:),dirCosinesx1(i,:));
end

% Element stiffness matrix

for k = 1:n % element stiffness matrix
    as = modElasticity(k)*area(k)/le(k);
    ts = shearModRig(k)*polMomInertia(k)/le(k);
    az1 = 12*modElasticity(k)*momInertiaZ(k)/le(k)^3;
    bz1 = 6*modElasticity(k)*momInertiaZ(k)/le(k)^3;
    cz1 = 4*modElasticity(k)*momInertiaZ(k)/le(k);
    dz1 = 2*modElasticity(k)*momInertiaZ(k)/le(k);
    ay1 = 12*modElasticity(k)*momInertiaY(k)/le(k)^3;
    by1 = 6*modElasticity(k)*momInertiaY(k)/le(k)^3;
    cy1 = 4*modElasticity(k)*momInertiaY(k)/le(k);
    dy1 = 2*modElasticity(k)*momInertiaY(k)/le(k);
    eleStiffDash(:,:,k) = [as 0 0 0 0 0 -as 0 0 0 0 0; % element stiffness matrix in local cordinates
                           0 az1 0 0 0 bz1 0 -az1 0 0 0 bz1;
                           0 0 ay1 0 -by1 0 0 0 -ay1 0 -by1 0;
                           0 0 0 ts 0 0 0 0 0 -ts 0 0;
                           0 0 -by1 0 cy1 0 0 0 by1 0 dy1 0;
                           0 bz1 0 0 0 cz1 0 -bz1 0 0 0 dz1;
                           -as 0 0 0 0 0 as 0 0 0 0 0;
                           0 -az1 0 0 0 -bz1 0 az1 0 0 0 -bz1;
                           0 0 -ay1 0 by1 0 0 0 cy1 0 -by1 0;
                           0 0 0 -ts 0 0 0 0 0 ts 0 0;
                           0 0 -by1 0 dy1 0 0 0 by1 0 cy1 0;
                           0 bz1 0 0 0 dz1 0 -bz1 0 0 0 cz1];

    transMat(:,:,k) = [dirCosinesx1(k,1) dirCosinesx1(k,2) dirCosinesx1(k,3) 0 0 0 0 0 0 0 0 0;
                   dirCosinesy1(k,1) dirCosinesy1(k,2) dirCosinesy1(k,3) 0 0 0 0 0 0 0 0 0;
                   dirCosinesz1(k,1) dirCosinesz1(k,2) dirCosinesz1(k,3) 0 0 0 0 0 0 0 0 0;
                   0 0 0 dirCosinesx1(k,1) dirCosinesx1(k,2) dirCosinesx1(k,3) 0 0 0 0 0 0;
                   0 0 0 dirCosinesy1(k,1) dirCosinesy1(k,2) dirCosinesy1(k,3) 0 0 0 0 0 0;
                   0 0 0 dirCosinesz1(k,1) dirCosinesz1(k,2) dirCosinesz1(k,3) 0 0 0 0 0 0;
                   0 0 0 0 0 0 dirCosinesx1(k,1) dirCosinesx1(k,2) dirCosinesx1(k,3) 0 0 0;
                   0 0 0 0 0 0 dirCosinesy1(k,1) dirCosinesy1(k,2) dirCosinesy1(k,3) 0 0 0;
                   0 0 0 0 0 0 dirCosinesz1(k,1) dirCosinesz1(k,2) dirCosinesz1(k,3) 0 0 0;
                   0 0 0 0 0 0 0 0 0 dirCosinesx1(k,1) dirCosinesx1(k,2) dirCosinesx1(k,3);
                   0 0 0 0 0 0 0 0 0 dirCosinesy1(k,1) dirCosinesy1(k,2) dirCosinesy1(k,3);
                   0 0 0 0 0 0 0 0 0 dirCosinesz1(k,1) dirCosinesz1(k,2) dirCosinesz1(k,3)];

     eleStiff(:,:,k) = transMat(:,:,k)'*eleStiffDash(:,:,k)*transMat(:,:,k); % element stiffness in global cordinates

end

% Global stiffness matrix asssembly
globStiff = zeros(max(max(DOF)),max(max(DOF)));
for i = 1:max(max(DOF))
    for j = 1:max(max(DOF))
        if i == j
            globStiff(i,j) = 0;
            for k = 1:size(DOFConnect,1)
                for l = 1:size(DOFConnect,2)
                    if DOFConnect(k,l) == i
                       globStiff(i,j) =  globStiff(i,j) + eleStiff(l,l,k);
                    end
                end
            end
        else
            globStiff(i,j) = 0;
          
            for k = 1:size(DOFConnect,1)
                 ldof1 = 0;
                 ldof2 = 0;  
                for l = 1:size(DOFConnect,2)
                        if DOFConnect(k,l) == i 
                            ldof1 = l;
                        elseif DOFConnect(k,l) == j
                            ldof2 = l;
                        end
                            if l == size(DOFConnect,2)
                                if  ldof1==0 || ldof2==0
                                 ldof1 = 0;
                                 ldof2 = 0;
                                else
                                    globStiff(i,j) =  globStiff(i,j) + eleStiff(ldof1,ldof2,k);
                                end
                            end
                    end
                end
            end
        end
end
finalStiffMat = globStiff;
disp('The Global Stiffness Matrix is');
disp(finalStiffMat);

% Element load vector
pointLoadMoment = zeros(max(max(DOF)),1); % Nodal loads
pointLoadMomentNode = [7 31];
pointLoadMomentValue = [3000 3000];
for k = 1:length(pointLoadMomentNode)
pointLoadMoment(pointLoadMomentNode(k)) = pointLoadMoment(pointLoadMomentNode(k)) + pointLoadMomentValue(k);
end

uniLoady1 = zeros(n,1); % Uniformly distributed load in element
uniLoady1Element = [1];
uniLoady1Value = [0];
for k = 1:length(uniLoady1Element)
uniLoady1(uniLoady1Element(k)) = uniLoady1(uniLoady1Element(k)) + uniLoady1Value(k);
end

uniLoadz1 = zeros(n,1); % Uniformly distributed load in element
uniLoadz1Element = [5 6 7 8];
uniLoadz1Value = [-20.833 27.778 20.833 -27.778];
for k = 1:length(uniLoadz1Element)
uniLoadz1(uniLoadz1Element(k)) = uniLoadz1(uniLoadz1Element(k)) + uniLoadz1Value(k);
end

for k = 1:n
eleForcedash(:,1,k) = [0;
                   uniLoady1(k)*le(k)/2;
                   uniLoadz1(k)*le(k)/2;
                   0;
                   -uniLoadz1(k)*le(k)^2/12;
                   uniLoady1(k)*le(k)^2/12;
                   0;
                   uniLoady1(k)*le(k)/2;
                   uniLoadz1(k)*le(k)/2;
                   0;
                   uniLoadz1(k)*le(k)^2/12;
                   -uniLoady1(k)*le(k)^2/12;];
eleForce(:,1,k) = transMat(:,:,k)'*eleForcedash(:,1,k); % convert force from local to global
end

globForce = zeros(max(max(DOF)),1);
for i = 1:max(max(DOF))
    for k = 1:size(DOFConnect,1)
        for l = 1:size(DOFConnect,2)
            if DOFConnect(k,l) == i
                globForce(i,1) = globForce(i,1) + eleForce(l,1,k);
            end
        end
    end
end

globForce = globForce + pointLoadMoment;
finalForceVec = globForce;
disp('The Global Force Vector is');
disp(finalForceVec);

%%% Constriants
% Constraint type b1Qp1 + b2Qp2 = bo
b1b2 = [1;1;1;1;1;1;
        1;1;1;1;1;1;
        1;1;1;1;1;1;
        1;1;1;1;1;1]; % coefficients of contraint (number of rows indicates a constraint equations)
bo = [0;0;0;0;0;0;
      0;0;0;0;0;0;
      0;0;0;0;0;0;
      0;0;0;0;0;0];
constDOF = [1;2;3;4;5;6;
            24;23;22;21;20;19;
            30;29;28;27;26;25;
            48;47;46;45;44;43
         ];
    ; % DOF constrained (one row one equation)
[row,cols] = size(constDOF);
C = max(max(globStiff))*10^8; % Penalty method


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

eqnn = length(globForce); % Number of equation

% Forward elimination (reduction of A,b)
for k = 1:eqnn-1
    for i = k+1:eqnn
        c = globStiff(i,k)/globStiff(k,k);
        for j=k:eqnn
           globStiff(i,j) = globStiff(i,j)-c*globStiff(k,j);
        end
        globForce(i) = globForce(i)-c*globForce(k);
    end
end


% Back substitution
globForce(eqnn) = globForce(eqnn)/globStiff(eqnn,eqnn);
for ii = 1:eqnn-1
    i = eqnn-ii;
    sum = 0;
    for j = i+1:eqnn ;
        sum = sum + globStiff(i,j)*globForce(j);
    end
    globForce(i) = (globForce(i)-sum)/globStiff(i,i);
end

for i=1:max(max(DOF))
finalDispVec(i,:) = [i globForce(i)];
end
disp('The displacements are')
disp(finalDispVec); % Vector giving displacement

% Shear Force And Bending Moment
for i = 1:n
    localDispVec(:,1,i) = transMat(:,:,i)*[finalDispVec(DOFConnect(i,1));
                                           finalDispVec(DOFConnect(i,2));
                                           finalDispVec(DOFConnect(i,3));
                                           finalDispVec(DOFConnect(i,4));
                                           finalDispVec(DOFConnect(i,5))
                                           finalDispVec(DOFConnect(i,6));
                                           finalDispVec(DOFConnect(i,7));
                                           finalDispVec(DOFConnect(i,8));
                                           finalDispVec(DOFConnect(i,9));
                                           finalDispVec(DOFConnect(i,10));
                                           finalDispVec(DOFConnect(i,11));
                                           finalDispVec(DOFConnect(i,12))];
end

for k = 1:n
eleEndEqLoad(:,:,k) = eleStiffDash(:,:,k) * localDispVec(:,1,k)  - eleForcedash(:,1,k);
end

for k = 1:n
% shear force and bending moment about local y axis
shearForcey1(k,:) = eleEndEqLoad(2,:,k); % Shear force at local node 1
shearForcey2(k,:) = -1*eleEndEqLoad(8,:,k); % Shear force at local node 2
bendingMomenty1(k,:) = -1*eleEndEqLoad(5,:,k); % Bending moment at local node 1
bendingMomenty2(k,:) = eleEndEqLoad(11,:,k); % bending moment at local node 2
shearForcey(k,:) = [shearForcey1(k,:) shearForcey2(k,:)]; % for an element
bendingMomenty(k,:) = [bendingMomenty1(k,:) bendingMomenty2(k,:)]; % for an element
end
disp('Shear force along local y axis is given by');
disp(shearForcey)
disp('Bending moment along local y axis is given by');
disp(bendingMomenty)

for k = 1:n
% shear force and bending moment about local z axis
shearForcez1(k,:) = eleEndEqLoad(3,:,k); % Shear force at local node 1
shearForcez2(k,:) = -1*eleEndEqLoad(9,:,k); % Shear force at local node 2
bendingMomentz1(k,:) = -1*eleEndEqLoad(6,:,k); % Bending moment at local node 1
bendingMomentz2(k,:) = eleEndEqLoad(12,:,k); % bending moment at local node 2
shearForcez(k,:) = [shearForcez1(k,:) shearForcez2(k,:)]; % for an element
bendingMomentz(k,:) = [bendingMomentz1(k,:) bendingMomentz2(k,:)]; % for an element
end
disp('Shear force along local z axis is given by');
disp(shearForcez)
disp('Bending moment along local z axis given by');
disp(bendingMomentz)


for k = 1:n
% Torsion about local x axis
torsionx1(k,:) = eleEndEqLoad(4,:,k); % torsion at local node 1
torsionx2(k,:) = eleEndEqLoad(10,:,k); % torsion at local node 2
torsionx(k,:) = [bendingMomentz1(k,:) bendingMomentz2(k,:)]; % for an element
end
disp('Torsion along local x axis given by');
disp(torsionx)

% Axial Force % -ve compression +ve tension
for i = 1:n
    AxialForce(i,1) = -eleEndEqLoad(1,:,i);
    Stress(i,1) = AxialForce(i,1)/area(i);
end
AxialForceStress = [AxialForce Stress];
disp('The axial force and stress in members')
disp(AxialForceStress)

% Reaction calculation R = KQ-F (only constrained DOF)
for i = 1:row
    for j = 1:cols
        k = constDOF(i,j);
reaction(i,1) = finalStiffMat(k,:)*finalDispVec(:,2) - finalForceVec(k);
    end
end
disp('The Reactions Are')
disp(reaction)

 





