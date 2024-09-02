clc
clear all
clear

% % % P-M Interaction curve ( Group 14 ) % % %

% Note : Computation time may be long please wait

% Section properties
D = 2200 ; % mm depth of section % % %
B = 2200; % mm breadth of section % % %
fck = 35 ; % MPa characteristic strength of concrete % % %
fy = 500 ; % MPa yield strength of steel % % %
rfdia = 16 ; % mm r/f bar diameter 
nRf = [10 10 4 4 4 4 14 14 4 4 4 4 10 10 ] ; % No of bars in each layer from top % % %
effCover = 56 ; % mm effective cover % % %
nSection = 11 ; % Divide s/c into discrete % % %
nSteel = 14 ; % No of layer of steel % % %
areaStrip = [280000 80000 80000 80000 80000 440000 80000 80000 80000 80000 280000] ; % Area of strip from top to bottom % % %

choice = menu ('Enter your choice', 'With partial factor of safety', 'Without partial factor of safety') ;
switch choice
    case 1 
        parSafconc = .45 ;
        parSafsteel = .87 ;
    case 2
        parSafconc = .67 ;
        parSafsteel = 1 ;
end



% 3 cases: xu >= D, D > xu >= 0, xu < 0 
% when s/c is in compression Emax = 0.0035 - 0.75*Emin

% xu >= D
k = 1;
o = 1;
l = 1;
for xu = 500000:-1:1
    neutralAxis(k) = xu ;
    if xu > D
        strainMat = [1 -xu/(xu-D); 1 .75] \ [0; 0.0035]; % Matrix containing [Emax; Emin]
        Emax = strainMat(1);
        Emin = strainMat(2);
        xCl = [100 300 500 700 900 1100 1300 1500 1700 1900 2100]; % Centreline of each strip from top % % %
        for i = 1:nSection
            strainCi(i) = ((xu - xCl(i)) / xu) * Emax; % Vector containing strain corresponding to centre line of each strip
        end

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

        lCi = (D / 2 - xCl) / 1000; % Vector containing distance in m of center line of strip from centroidal axis

        % Vector containing moment due to each strip
        mCi = pCi .* lCi;

        TotalMconc = sum(mCi); % Total moment on concrete section (kN*m)


        % Steel capacity
        xSl(1) = effCover; % Centreline distance of bar from top % % %
        xSl(length(nRf)) = D - effCover ;
        distbtwbar = (D - effCover * 2 ) / (length(nRf)-1) ; % distance between two layers of steel
        for i = 1:length(nRf)-1
            xSl (i+1) = xSl(i) + distbtwbar ;
        end


        for i = 1:nSteel
            strainSi(i) = ((xu - xSl(i)) / xu) * Emax; % Vector containing strain corresponding to r/f
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

        % Vector containing distance in m of center line of steel from centroidal axis
        lSi = (D / 2 - xSl) / 1000;

        % Vector containing moment due to each steel layer
        mSi = pSi .* lSi;

        TotalMSteelc = sum(mSi); % Total moment due to steel (kN*m)

        totalLoadCap = totalPConc + totalPSteel; % Total load capacity of s/c (kN)
        totalMomCap = TotalMconc + TotalMSteelc; % Total moment capacity of s/c (kN*m)

        P1(o) = totalLoadCap;
        M1(o) = totalMomCap;

        o = o + 1;




    elseif xu == D
        Emax = 0.0035;
        Emin = 0;

        xCl = [100 300 500 700 900 1100 1300 1500 1700 1900 2100] ; % Vector containing centre line distance of each strip from top of section % % %

        for i = 1:length(xCl)
            strainCi(i) = ((xu - xCl(i)) / xu) * Emax; % Vector containing strain corresponding to centre line of each strip
        end

        for i = 1:length(strainCi) % Vector containing stress corresponding to strain at centre line of strip
            if strainCi(i) < 0.002
                stressCi(i) = parSafconc * fck * (2 * (strainCi(i) / 0.002) - (strainCi(i) / 0.002)^2) ;
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

        lCi = (D / 2 - xCl) / 1000;

        % Vector containing moment due to each strip
        mCi = pCi .* lCi;

        TotalMconc = sum(mCi); % Total moment on concrete section (kN*m)

        % Steel capacity
        xSl(1) = effCover; % Centreline distance of bar from top % % %
        xSl(length(nRf)) = D - effCover ;
        distbtwbar = (D - effCover * 2 ) / (length(nRf)-1) ; % distance between two layers of steel
        for i = 1:length(nRf)-1
            xSl (i+1) = xSl(i) + distbtwbar ;
        end



        for i = 1:nSteel
            strainSi(i) = ((xu - xSl(i)) / xu) * Emax; % Vector containing strain corresponding to r/f
        end

        E1 = 0.0018; % Nonlinear part of stress-strain curve
        E2 = ((parSafsteel * fy) / (2 * 10^5)) + 0.002;
        S1 = 350;
        S2 = parSafsteel * fy;
        for i = 1:nSteel % Vector containing stress in steel corresponding to strain
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
        end

        % Vector containing force in each strip in kN
        pSi = stressSi .* ((pi / 4) * rfdia^2 .* nRf / 1000);

        totalPSteel = sum(pSi); % Total compression force on steel (kN)

        % Vector containing distance in m of center line of steel from centroidal axis
        lSi = (D / 2 - xSl) / 1000;

        % Vector containing moment due to each steel layer
        mSi = pSi .* lSi;

        TotalMSteelc = sum(mSi); % Total moment due to steel (kN*m)

        totalLoadCap = totalPConc + totalPSteel; % Total load capacity of s/c (kN)
        totalMomCap = TotalMconc + TotalMSteelc; % Total moment capacity of s/c (kN*m)

        P3 = totalLoadCap;
        M3 = totalMomCap;
    clear strainCi stressCi strainSi stressSi xCl xSl pCi mCi pSi mSi


    else % xu < 0
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
        for i = 1:length(xCl)
            strainCi(i) = ((xu-xCl(i))/xu)*0.0035; % Strain corresponding to the centerline of each strip
        end

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

            lCi = (D/2 - xCl)/1000; % Distance of the centerline of the strip from the centroidal axis (m)

            mCi = pCi.*lCi ; % Moment due to each strip (kNm)

            TotalMconc = sum(mCi) ; % Total moment on the concrete section (kNm)
        else
            Beff = 1400 ;  % xu at top strip
            pCi(1) =(stressCi(1)*xu*Beff)/1000;
            mCi(1) = pCi(1)*(D/2-xCl(1))/1000;
            totalPConc = pCi(1);
            TotalMconc = mCi(1); % Total moment on the concrete section (kNm)
        end

        % Steel capacity
        xSl(1) = effCover; % Centreline distance of bar from top % % %
        xSl(length(nRf)) = D - effCover ;
        distbtwbar = (D - effCover * 2 ) / (length(nRf)-1) ; % distance between two layers of steel
        for i = 1:length(nRf)-1
            xSl (i+1) = xSl(i) + distbtwbar ;
        end

        strainSi = ((xu - xSl) ./ xu) .* 0.0035; % Vector containing strain corresponding to r/f



        E1 = 0.0018; % nonlinear part of stress-strain curve
        E2 = ((parSafsteel * fy) / (2 * 10^5)) + 0.002;
        S1 = 350 ;
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

        % Vector containing distance in m of center line of steel from centroidal axis
        lSi = (D/2 - xSl) / 1000;

        % Vector containing moment due to each steel layer
        mSi = pSi .* lSi;

        TotalMSteelc = sum(mSi); % Total moment due to steel in kNm

        totalLoadCap = totalPConc + totalPSteel; % Total load capacity of s/c
        totalMomCap = TotalMconc + TotalMSteelc;  % Total moment capacity of s/c

        P2(l) = totalLoadCap;
        M2(l) = totalMomCap;
        l = l + 1;
        clear strainCi stressCi strainSi stressSi xCl xSl pCi mCi pSi mSi lCi lSi

    end

    k = k  + 1 ;

  end

% Final data for plotting PM curve
mFinal = [M1 M3 M2] ;
pFinal = [P1 P3 P2] ;



% Plot excel point on the graph (Verification)
format longg
xuValue1 = 4400 ;
xuValue2 = 2200 ;
xuValue3 = 1100 ; 
verifyXu = [xuValue1, xuValue2, xuValue3] ;
indexV1 = find (neutralAxis == xuValue1) ;
indexV2 = find (neutralAxis == xuValue2) ;
indexV3 = find (neutralAxis == xuValue3) ;
point1 = [mFinal(indexV1) pFinal(indexV1)] ;  % xu = 4400 mm
point2 = [mFinal(indexV2) pFinal(indexV2)] ;  % xu = 2200 mm
point3 = [mFinal(indexV3) pFinal(indexV3)] ;  % xu = 1100 mm
verifyM = [point1(1) point2(1) point3(1)] ;
verifyP = [point1(2) point2(2) point3(2)] ;
color = {'r' ,'g', 'b'} ;

for i = 1:length(verifyP)
plot (verifyM(i) , verifyP(i), '*', Color=color{i}) ;
end



for i = 1:length(verifyP)
    text ( verifyM(i), verifyP(i), ['(', num2str(verifyM(i)), ',', num2str(verifyP(i)), ')'], 'HorizontalAlignment', 'left');
end
hold on 
     

% Plot PM interation curve
plot (mFinal, pFinal, 'r') ;
title ('P-M Interaction curve') ;
xlabel (' Moment (KNm) ')
ylabel (' Axial force (KN) ');
grid on
axis tight

   
            
