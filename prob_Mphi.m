clear all
clc

% % % M - Phi interaction curve ( Group 14 ) % % %

% Note : Computation time may be long please wait

% axialP = 3500 ; % KN Axial force 3500 KN and 0 KN
fck = 35 ; % MPa characteristic strength of concrete % % %
fy = 500 ; % MPa yield strength of steel % % %
B = 2200 ; % width of section % % %
D = 2200;   % D of section % % %
areaEntire = (B * D) - (4 * 1000 * 800) ; % Area of entire section % % %
% phi = 1e-5; % per m assume the curvature
effCover = 56 ; % Effective cover % % %
nSteel = 14 ; % No of layer of steel % % %
rfdia = 16; % mm r/f bar diameter
nRf = [10 10 4 4 4 4 14 14 4 4 4 4 10 10 ] ; % No of bars in each layer from top % % %
nSection = 11 ; % Divide s/c into discrete % % %
areaStrip = [280000 80000 80000 80000 80000 440000 80000 80000 80000 80000 280000] ; % Area of strip from top to bottom % % %


% Initial calculation

choice2 = menu ('Enter your choice', 'With partial factor of safety', 'Without partial factor of safety') ;
switch choice2
    case 1 
        parSafconc = .45 ;
        parSafsteel = .87 ;
    case 2
        parSafconc = .67 ;
        parSafsteel = 1 ;
end

choice1 = menu('Select required axial force', '0 KN', '3500 KN') ;
switch choice1
    case 1
        axialP = 0;
        finalStrain = .0002 ;
        totalP = 100 ; % Initialize totalP to enter into the loop
    case 2
        axialP = 3500;
        stressC = axialP * 1000 / areaEntire ; % Find the uniform stress in the section
        coefficientPoly = [-(1/.002)^2 * parSafconc * fck, parSafconc * fck * 2 / .002, -stressC] ;
        rootsPoly = roots (coefficientPoly) ; % strain corresponding to uniform stress
        initialStrain = min( abs(rootsPoly)) ;
        finalStrain = initialStrain ; % Initialize final strain
        totalP = 0 ; % Initialize totalP to enter into the loop
end



o = 1; 

for phi = 4e-5 : 1e-4 : .02 % loop for different curvature
while totalP  > axialP+20 || totalP < axialP-20
    clear strainCi stressCi strainSi stressSi xCl xSl pCi mCi pSi mSi
if totalP  > axialP+20
    finalStrain = finalStrain * (2/3) ; % Dncrease the strain since Pcap < axialP
else
    finalStrain = finalStrain * 3.4722 ; % Iecrease the strain since Pcap < axialP
end
xu = (finalStrain / phi) *1000 ;  % Neutral axis corresponding to finalStrain

if xu >= D
         
        xCl = [100 300 500 700 900 1100 1300 1500 1700 1900 2100]; % Centreline of each strip from top % % %
        strainCi = (( xu - xCl ) .* finalStrain) ./ xu ; % Strain in centreline of strip
        for i = 1:length(strainCi) % Vector containing stress corresponding to strain at centre line of strip
            if strainCi(i) < 0.002
                stressCi(i) = parSafconc * fck * (2 * (strainCi(i) / 0.002) - (strainCi(i) / 0.002)^2);
            elseif strainCi(i) >= 0.002
                stressCi(i) = parSafconc * fck;
            else
                fprintf('Error');
            end
        end

        % Concrete capacity
        % Vector containing force in each strip in kN % % %
        pCi = stressCi .* areaStrip / 1000; 

        totalPConc = sum(pCi); % Total compression force on concrete section (kN)

        % Steel capacity
        xSl(1) = effCover; % Centreline distance of bar from top % % %
        xSl(length(nRf)) = D - effCover ;
        distbtwbar = (D - effCover * 2 ) / (length(nRf)-1) ; % distance between two layers of steel
        for i = 1:length(nRf)-1
            xSl (i+1) = xSl(i) + distbtwbar ;
        end

        for i = 1:nSteel
            strainSi(i) = ((xu - xSl(i)) / xu) * finalStrain ; % Vector containing strain corresponding to r/f
        end

        E1 = 0.0018; % Nonlinear part of stress-strain curve
        E2 = ((parSafsteel * fy) / (2 * 10^5)) + 0.002;
        S1 = 350;
        S2 = parSafsteel * fy;
        for i = 1:length(strainSi) % Vector containing stress in steel corresponding to strain
            if strainSi(i) <= E1
                stressSi(i) = S1 * strainSi(i) / E1;
            elseif strainSi(i) > E1 && strainSi(i) <= E2
                a = (S2 - S1) / ((E2^2 - E1^2) - (2 * E2 * (E2 - E1)));
                b = -2 * a * E2;
                c = S1 - a * E1^2 - b * E1;
                stressSi(i) = a * strainSi(i)^2 + b * strainSi(i) + c;
            else
                stressSi(i) = parSafsteel * fy;
            end
        end

        % Vector containing force in each strip in kN
        pSi = stressSi .* ((pi / 4) * rfdia^2 .* nRf / 1000);

        totalPSteel = sum(pSi); % Total compression force on steel (kN)

        totalP = totalPConc + totalPSteel ;

elseif xu < D
        % Calculate the centerline distance of each strip from the top of the section
      % Calculate the centerline distance of each strip from the top of the section
        Strip1 = (D/nSection);
        if Strip1 < xu
            xCl(1) = (D/nSection)/2;
            for i = 1:1000
                if xCl(i)+(D/nSection)*(3/2) <= xu
                    xCl(i+1) = xCl(i)+D/nSection;
                elseif xCl(i)+(D/nSection)/2 == xu
                    break ;
                else 
                    xCl(i+1) = xCl(i)+(D/nSection)/2+(xu-(xCl(i)+(D/nSection)/2))/2;
                    break;
                end
            end
        else 
            xCl(1) = xu/2;
        end
        
        % Calculate strain for each strip
        strainCi = (( xu - xCl ) .* finalStrain) ./ xu ; 
        
        % Calculate stress for each strip
        for i = 1:length(xCl)
            if strainCi(i) < 0.002
                stressCi(i) =  parSafconc * fck * (2 * (strainCi(i) / 0.002) - (strainCi(i) / 0.002)^2);
            elseif strainCi(i) >= 0.002
                stressCi(i) = parSafconc*fck;
            else 
                fprintf('Error');
            end
        end
        
        % Concrete capacity
        if Strip1 < xu
            for i = 1:length(xCl)

            if xCl(i) <= 2200 && xCl(i) > 2000
                Beff = 1400 ;
            elseif xCl(i) <= 2000 && xCl(i) > 1200
                Beff = 400 ;
            elseif xCl(i) <= 1200 && xCl(i) > 1000
                Beff = 2200 ;   
            elseif xCl(i) <= 1000 && xCl(i) > 200
                Beff = 400 ;
            elseif xCl(i) <= 200 && xCl(i) >= 0
                Beff = 1400 ;
            else
                printf ( 'xu < 0') ;
            end
            if i == length (xCl)
                break ;
            end
                pCi(i) = (stressCi(i)*(Beff*(D/nSection)))/1000; % Force in each strip in kN considering effective width ( reduced width ) % % % 
            end 
            pCi(i) = (stressCi(i)*(xu-(xCl(i-1)+Strip1/2))*Beff)/1000;

            totalPConc = sum(pCi); % Total compression force on the 2200 x 2200 concrete section (kN)

       else
            Beff = 1400 ;  % xu at top strip
            pCi(1) =(stressCi(1)*xu*Beff)/1000;
            totalPConc = pCi(1);
        end

        
         
        % Steel capacity
        xSl(1) = effCover; % Centreline distance of bar from top % % %
        xSl(length(nRf)) = D - effCover ;
        distbtwbar = (D - effCover * 2 ) / (length(nRf)-1) ; % distance between two layers of steel
        for i = 1:length(nRf)-1
            xSl (i+1) = xSl(i) + distbtwbar ;
        end
        
        strainSi = ((xu - xSl) ./ xu) .* finalStrain; % Vector containing strain corresponding to r/f
           
               
        
        E1 = 0.0018; % nonlinear part of stress-strain curve
        E2 = ((parSafsteel * fy) / (2 * 10^5)) + 0.002;
        S1 = 350;
        S2 = parSafsteel * fy;
        
        for i = 1:nSteel  % Vector containing stress in steel corresponding to strain
            if strainSi(i) > 0
                if strainSi(i) < E1
                    stressSi(i) = S1 * strainSi(i) / E1;
                elseif strainSi(i) > E1 && strainSi(i) < E2
                    a = (S2 - S1) / ((E2^2 - E1^2) - (2 * E2 * (E2 - E1)));
                    b = -2 * a * E2;
                    c = S1 - a * E1^2 - b * E1;
                    stressSi(i) = a * strainSi(i)^2 + b * strainSi(i) + c;
                else
                    stressSi(i) = parSafsteel * fy;
                end
            elseif strainSi(i) < 0
                if abs(strainSi(i)) < E1
                    stressSi(i) = -(S1 * abs(strainSi(i)) / E1);
                elseif abs(strainSi(i)) > E1 && abs(strainSi(i)) < E2
                    a = (S2 - S1) / ((E2^2 - E1^2) - (2 * E2 * (E2 - E1)));
                    b = -2 * a * E2;
                    c = S1 - a * E1^2 - b * E1;
                    stressSi(i) = -(a * abs(strainSi(i))^2 + b * abs(strainSi(i)) + c);
                else
                    stressSi(i) = -parSafsteel * fy;
                end
            else
                stressSi(i) = 0;
            end
        end
        
        % Vector containing force in each strip in kN
        pSi = (stressSi .* ((pi/4) * rfdia^2 .* nRf) / 1000);
       
        
        totalPSteel = sum(pSi); % Total compression force on steel in kN
        totalP = totalPConc + totalPSteel ;

       
else
        printf (" xu < 0 ") ;
end

    
end




% Calculate moment capacity of section
        
        neutralAxis (o) = xu;

        % Concrete capacity
        lCi = (D / 2 - xCl) / 1000; % Vector containing distance in m of center line of strip from centroidal axis

        % Vector containing moment due to each strip
        mCi = pCi .* lCi;

        TotalMconc = sum(mCi); % Total moment on concrete section (kN*m)
         
        % Steel capacity
        % Vector containing distance in m of center line of steel from centroidal axis
        lSi = (D / 2 - xSl) / 1000;

        % Vector containing moment due to each steel layer
        mSi = pSi .* lSi;

        TotalMSteelc = sum(mSi); % Total moment due to steel (kN*m)

        TotalM = TotalMconc + TotalMSteelc ;
        

        curv(o) = phi ;
        Mn(o) = TotalM;

        o = o + 1 ;

        if axialP == 0
           totalP = 100 ;
        else
           totalP = 0 ;
        end

end

% plot M-Phi curve
plot ( curv, Mn, 'r') ;
title ('M - Phi curve') ;
xlabel ('curvature (/m)') ;
ylabel (' Moment (KNm) ') ;
grid on ;
axis tight ;









