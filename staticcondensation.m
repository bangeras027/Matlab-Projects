function [freq,reArrModeShape] = staticcondensation(stiffMat, massMat, primaryDof, secondaryDof)
% Static condensation method to reduced DOF in dynamic problemstn
% Input the full stiffness and full mass matrix, outputs the natural
% frequency and mode shape including the condensed dof
% stiffMat = stiffness matrix to condense
% massMat = mass matrix to condense
% primary = row vector containing no condensing dof
% secondary = row vector containing dof to be condensed
% output is natural frequency and modeshape

dofOrder = [secondaryDof primaryDof] ; % DOF order secondary dof first followed by primary dof

% Condensation process
for i = 1:size(stiffMat,1) % Modified mass matrix
    for j = 1:size(stiffMat,1)
        modMassMat(i,j) = massMat(dofOrder(i),dofOrder(j)) ;
    end
end

for i = 1:size(stiffMat,1) % Modified stiffness matrix
    for j = 1:size(stiffMat,1)
        modStiffMat(i,j) = stiffMat(dofOrder(i),dofOrder(j)) ;
    end
end

for i = 1:size(stiffMat,1) % Gauss elimination method to get the transformation matrix and modified stiff matrix
    if i <= length (secondaryDof) 
        modStiffMat(i,:) = modStiffMat(i,:) ./ modStiffMat(i,i) ;
        for k = i+1:size(stiffMat,1)
        modStiffMat(k,:) = modStiffMat(k,:) - modStiffMat(k,i) * modStiffMat(i,:) ;
        end
    end
end  

Tdash = -1*modStiffMat(1:length(secondaryDof),length(secondaryDof)+1:size(stiffMat,1)) ; 
condensedStiff = modStiffMat(length(secondaryDof)+1:size(stiffMat,1),length(secondaryDof)+1:size(stiffMat,1)) ; % Condensed stiffness matrix
T = [ Tdash; eye(size(Tdash,2))] ; % Transformation matrix
condensedMass = T' * massMat * T ; % Condensed mass matrix
[eigenVector, eigenValue] = eig(condensedStiff,condensedMass) ; % Eigen value problem
omega = sqrt(diag(eigenValue)) ; % Circular frequency (rad/s)
freq = omega/ (2*pi) ; % Natural frequency (cps)
for i = 1:size(eigenVector, 2) % Normalize the eigenvectors with respect to the mass matrix
    eigenVector(:, i) = -1*eigenVector(:, i) / sqrt(eigenVector(:, i)' * condensedMass * eigenVector(:, i));
end

for i = 1:size(stiffMat,1)-length(secondaryDof)
modeShape(:,i) = T * eigenVector(:,i) ;
end

for i = 1:size(stiffMat,1) % Rearrange the mode shapes
    hold = dofOrder(i) ;
    reArrModeShape(hold,:) = modeShape(i,:) ;
end
end