function [freq,reArrModeShape] = dynamiccondensation(stiffMat, massMat, primaryDof, secondaryDof)
% Dynamic condensation method to reduced DOF in dynamic problemstn
% Input the full stiffness and full mass matrix, outputs the natural
% frequency and mode shape including the condensed dof
% stiffMat = stiffness matrix to condense
% massMat = mass matrix to condense
% primary = row vector containing no condensing dof
% secondary = row vector containing dof to be condensed
% output is natural frequency and modeshape

dofOrder = [secondaryDof primaryDof] ; % DOF order secondary dof first followed by primary dof

% Condensation process
for i = 1:size(stiffMat,2) % Modified mass matrix
    for j = 1:size(stiffMat,2)
        modMassMat(i,j) = massMat(dofOrder(i),dofOrder(j)) ;
    end
end

for i = 1:size(stiffMat,2) % Modified stiffness matrix
    for j = 1:size(stiffMat,2)
        modStiffMat(i,j) = stiffMat(dofOrder(i),dofOrder(j)) ;
    end
end

omegaSq(1) = 0 ; % Initialize trial value for first fundamental natural frequency
tol = 100 ; % Enther the while loop

for l = 1:size(stiffMat,2)-length(secondaryDof)
    while tol > .01
    dynamicMat = modStiffMat - omegaSq(l) * modMassMat ;

    for i = 1:size(stiffMat,2) % Gauss elimination method to get the transformation matrix and modified stiff matrix
        if i <= length (secondaryDof) 
            dynamicMat(i,:) = dynamicMat(i,:) ./ dynamicMat(i,i) ;
            for k = i+1:size(stiffMat,2)
            dynamicMat(k,:) = dynamicMat(k,:) - dynamicMat(k,i) * dynamicMat(i,:) ;
            end
        end
    end 
    
    Tdash = -1*dynamicMat(1:length(secondaryDof),length(secondaryDof)+1:size(stiffMat,2)) ; 
    condensedDynamicMat = dynamicMat(length(secondaryDof)+1:size(stiffMat,2),length(secondaryDof)+1:size(stiffMat,2)) ; % Condensed stiffness matrix
    T = [ Tdash; eye(size(Tdash,2))] ; % Transformation matrix
    condensedMass = T' * massMat * T ; % Condensed mass matrix
    condensedStiff = condensedDynamicMat + omegaSq(l) * condensedMass ; % Condensed stiff matrix
    [eigenVector, eigenValue] = eig(condensedStiff,condensedMass) ; % Eigen value problem
    for i = 1:size(eigenVector, 2)
        eigenVector(:, i) = -1*eigenVector(:, i) / sqrt(eigenVector(:, i)' * condensedMass * eigenVector(:, i));
    end
    [eigenValue, sortIndex] = sort(diag(eigenValue)) ; % Circular frequency (rad/s)
    eigenVector = eigenVector(:, sortIndex) ;
    tol = abs(eigenValue(l)-omegaSq(l)) ;

    if l < length(primaryDof)
        omegaSq(l) = eigenValue(l) ;
        omegaSq(l+1) = eigenValue(l+1) ;
    else
        omegaSq(l) = eigenValue(l) ;
    end
    
    eigenVector = eigenVector(:,l) ;
    end
    
    omega(l) = sqrt(omegaSq(l)) ; 
    tol = 100 ;
    eigenVector = -eigenVector(:,1) ;
    modeShape = T * eigenVector ;
    
    for i = 1:size(stiffMat,2) % Rearrange the mode shapes
        hold = dofOrder(i) ;
        reArrModeShape(hold,l) = modeShape(i) ;
    end
   freq = (omega)/ (2*pi) ; 
end
end

