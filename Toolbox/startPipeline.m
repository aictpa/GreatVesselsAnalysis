function out = startPipeline(datasetPath,current_path,showsFigures)
% STARTPIPELINE  Starting model pipeline
%
%   Examples:
%      out = STARTPIPELINE(datasetPath,current_path,showsFigures)
%      "datasetPath" holds location of the data.
%      If the "showsFigures" is set to "true", the output of the task will
%      be illustrated.

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


out = false;

f = dir(datasetPath);
exams_list = fullfile(datasetPath,{f(3:end).name});
iScans = size(exams_list,2);


%% Step 1 setting up variables

imgTemp = cell(1,iScans);
lastDesAorta = cell(1,iScans);
lastDesAortaLoc = cell(1,iScans);
ascAorta = cell(1,iScans);
ascLoc = cell(1,iScans);
imgTraB = cell(1,iScans);
imgTraO = cell(1,iScans);
imgBifB = cell(1,iScans);
imgBifO = cell(1,iScans);
imgCarinaB = cell(1,iScans);
imgCarinaO = cell(1,iScans);
imgDesaC = cell(1,iScans);
imgDesaB = cell(1,iScans);
imgDesaO = cell(1,iScans);
imgAscaB = cell(1,iScans);
imgAscaO = cell(1,iScans);
imgPTB = cell(1,iScans);
imgPTB2 = cell(1,iScans);
imgPTO = cell(1,iScans);
imgPVO = cell(1,iScans);
imgPTDia = cell(1,iScans);
imgPTHUMask = cell(1,iScans);
imgPTHU = cell(1,iScans);
representationOfApicalLevelPVLoc = cell(1,iScans);
representationOfApicalLevelPV = cell(1,iScans);
seedPointOfTra = cell(1,iScans);
seedPointOfBifR = cell(1,iScans);
seedPointOfBifL = cell(1,iScans);
seedPointOfBif = cell(1,iScans);
seedPointOfDesA = cell(1,iScans);
seedPointOfAscA = cell(1,iScans);
seedPointOfPT = cell(1,iScans);
seedPointOfPV = cell(1,iScans);
idxOfPre = cell(1,iScans);
indexRho = cell(1,iScans);
distanceCm = cell(1,iScans);
rightMidPoint = cell(1,iScans);
leftMidPoint = cell(1,iScans);
ptCan = cell(1,iScans);
ptLoc = cell(1,iScans);
ptCir = cell(1,iScans);
ptCirD = cell(1,iScans);
ptDia = cell(1,iScans);
ptDia3 = cell(1,iScans);
ptDiaLines  = cell(1,iScans);
ptDiaLines3  = cell(1,iScans);
ptSize = cell(1,iScans);
hu = cell(1,iScans);
huN = cell(1,iScans);
ascDia = cell(1,iScans);
traDia = cell(1,iScans);
rightCentroids = cell(1,iScans);
leftCentroids = cell(1,iScans);
l12c = cell(1,iScans);
SDdesA = cell(1,iScans);
desAorta = cell(1,iScans);
SDdesALoc = cell(1,iScans);
rotated1=cell(1,iScans);
rotated2=cell(1,iScans);
rotated3=cell(1,iScans);
ftimeMain = cell(1,iScans);
ftimeImgUpload = cell(1,iScans);
ftimeOriControl = cell(1,iScans);
ftimeCarina  = cell(1,iScans);
ftimeDesA = cell(1,iScans);
ftimeAscA = cell(1,iScans);
ftimePT = cell(1,iScans);

counter = 0;
colorMap=[1 0 0; 0 0 1];

ts=cputime;
tic

for iCase = 1:iScans
    try        
               
        tStartMain = tic;     
        
        file_list = regexp(exams_list(iCase), filesep, 'split');
        timestamp_str = datestr(now,'yyyy_mm_dd_HH_MM_SS_FFF');
        
        %% Step 1 Converting CT Data to Hounsfield Units
        tStartImgUpload = tic;
        [~,sliceHU,sliceInfo] = getCTPAExam(exams_list{iCase});
        ctpa = reshape([sliceHU(1:end).HU],[size(sliceHU(1).HU),length(sliceHU)]);
        ftimeImgUpload{iCase} = toc(tStartImgUpload);
        
        
        %% Step 2 Finding scanning direction
        tStartOriControl= tic;
        sizeOfRawCTPA = size(ctpa,3);
        if ( sliceInfo(2).info.ImagePositionPatient(3) > sliceInfo((end-1)).info.ImagePositionPatient(3) )
            PATIENT_POSITION = 1; % first starts with top slice            
            rotationControlIndex = 10;
            reverseRotationControlIndex = sizeOfRawCTPA-1;
            imgIN =  ctpa(:,:,rotationControlIndex); % 2
            loopIncrementVal = 1;
            startIndex = 1;            
        else % idxOfTop == 2 starts with bottom slice
            PATIENT_POSITION = 2;            
            rotationControlIndex = size(ctpa,3) - 10;
            reverseRotationControlIndex = 2;
            imgIN =  ctpa(:,:,rotationControlIndex); % -1
            loopIncrementVal = -1;
            startIndex = sizeOfRawCTPA;            
        end
                
        pixelSpacing = sliceInfo(20).info.PixelSpacing(1);
        
        
        %% Step 3 Check Orientation
        
        rotationAngle = abs(calculateExamOrientation(imgIN));
        
        isRotated = false;
        
        % Rotate study
        if(rotationAngle ~= 0)
            
            %% FIX
            angle = findAccurateAngle(reverseRotationControlIndex,ctpa);
            
            %% Take care of zero pixel
            t = false(size(ctpa));
            t(ctpa == 0) = 1;
            t = imrotate(t,angle,'nearest','crop');
            
            ctpa = imrotate(ctpa,angle,'nearest','crop');
            
            ctpa(ctpa==0) = -1024;
            ctpa(t) = 0;
            
            isRotated = true;
            rotated1{iCase} = iCase;
            
        end
        
        ftimeOriControl{iCase} = toc(tStartOriControl);
          
        
        %% Step 4 Trachea Detection
        
        tStartCarina= tic;
        
        [tracheaS,sliceNoOfTrachea,seedPointOfTrachea,~] = tracheaDetection(ctpa,startIndex,loopIncrementVal);
      
       
        %% Step 5 Checking Tracheal Intubation
        
        [tracheaS,sliceNoOfTrachea,seedPointOfTrachea] = checkTrachealIntubation(ctpa,tracheaS,sliceNoOfTrachea,seedPointOfTrachea,PATIENT_POSITION,iCase);
                
        seedPointOfTra{iCase} = seedPointOfTrachea;        
        
        imgTraB{iCase} = tracheaS;
        imgTraO{iCase} = ctpa(:,:,sliceNoOfTrachea);
        
        % Diameter of the Trachea (Transverse/Lateral)
        stats_trachea = regionprops(tracheaS,'EquivDiameter');
        equivDia = stats_trachea(1).EquivDiameter;
        traDia{iCase} = equivDia*pixelSpacing;
        
        %% Step 6 Airways Segmentation        
        
        [imgTemp{iCase},bifurPoint,rightMainBronchus,leftMainBronchus,~,~,~,tresholdVal] = airwaysExtraction(ctpa,PATIENT_POSITION,seedPointOfTrachea,sliceNoOfTrachea);
        
        rightMainBronchus = bwareafilt(rightMainBronchus,1,8);
        leftMainBronchus = bwareafilt(leftMainBronchus,1,8);
        
        rightMainTracheaSeedPoints = getSeedPoint(imfill(rightMainBronchus,'holes'));
        
        leftMainTracheaSeedPoints = getSeedPoint(imfill(leftMainBronchus,'holes'));
        
        rightMainBronchusBif = rightMainBronchus;
        leftMainBronchusBif = leftMainBronchus;
        
        imgBifB{iCase} = rightMainBronchusBif;%
        imgBifO{iCase} = ctpa(:,:,bifurPoint);%
        
        
        seedPointOfBifR{iCase} = round(rightMainTracheaSeedPoints);
        seedPointOfBifL{iCase} = round(leftMainTracheaSeedPoints);
        seedPointOfBif{iCase} = [seedPointOfBifR{iCase}; seedPointOfBifL{iCase}];              
        
        
        %% Step 7 Find Possible Carina Level 
        
        [imgTemp{iCase},rightCentroids{iCase},leftCentroids{iCase},rightMainBronchus,leftMainBronchus,sliceLoc,~,rightMidPoint{iCase},leftMidPoint{iCase},distanceCm{iCase}] = carinaLevelDetection(ctpa,PATIENT_POSITION,pixelSpacing,rightMainBronchus,leftMainBronchus,bifurPoint,tresholdVal);
               
        
        tmp = imadd(rightMainBronchus,leftMainBronchus);
        tmp = bwlabeln(logical(tmp));          
        
        rightMainBronchusBeforeOrientation = rightMainBronchus;
                
        imgCarinaB{iCase} = rightMainBronchusBeforeOrientation;%
        imgCarinaO{iCase} = ctpa(:,:,sliceLoc);%
        
        angle = 0;
        addedAngle = 0;
        
        if (~isRotated)
            [imgTemp{iCase},angle] = singleRotation(imgTemp{iCase});
            
            if(abs(angle) > 25)                               
                
                leftMainBronchus = imrotate(leftMainBronchus,angle,'nearest','crop');
                leftCentroids{iCase} = getSeedPoint(imfill(leftMainBronchus,'holes'));
                
                rightMainBronchus = imrotate(rightMainBronchus,angle,'nearest','crop');
                rightCentroids{iCase} = getSeedPoint(imfill(rightMainBronchus,'holes'));
                                
                rotated2{iCase} = iCase;
                
                isRotated = true;
            end
            
            if (~isRotated)
                stats = regionprops(leftMainBronchus,'Orientation');               
                
                if(abs(stats.Orientation)>28 && abs(stats.Orientation)<45)
                    rotated3{iCase} = iCase;
                    %isRotated = true;
                    addedAngle = abs(stats.Orientation)-5;                    
                end
            end
            
        end % end of if
        
        
        ftimeCarina{iCase} = toc(tStartCarina);        
        
        
        %% Step 8 Detect Descending Aorta
        
        tStartDesA = tic;        
        
        SIGMA = 1.5;
        I = imgaussfilt(imgTemp{iCase},SIGMA);
        L = eigenvalHessian(I);
        
        E = edge(I,'canny');
        E = imdilate(E,strel('Disk',1));
        LE = (L + E) > 0;
        intensityOfHessian = sum(LE(:));
        l12c{iCase} = angle;%sum(L(:));
        if(intensityOfHessian < 49000)
            LE = imclose(LE,strel('disk',2));
        end
              
        [imgDesAorta,~,~,~]  = desAortaExtraction(imgTemp{iCase},I,LE,pixelSpacing,rightCentroids{iCase},leftCentroids{iCase}, addedAngle);
        
        
        %% Step 9 Calculate Density SD of Descending Aorta
        
        SDdesALoc{iCase} = sliceLoc;
        [SDdesA{iCase},desAorta{iCase},BWt,seedPoint,deviation] = getDesAortaSD(imgDesAorta,angle,ctpa,sliceLoc,pixelSpacing);
        
        imgDesaB{iCase} = desAorta{iCase};%
        imgDesaO{iCase} = ctpa(:,:,sliceLoc);%
        seedPointOfDesA{iCase} = seedPoint;
                
        ftimeDesA{iCase} = toc(tStartDesA);
        
        
        %% Step 10 Detect Aortic Arch
        tStartAscA = tic;
        
        [lastDesAorta{iCase},lastDesAortaLoc{iCase},BW,~] = aorticarchDetection(sliceLoc,desAorta{iCase},loopIncrementVal,ctpa,pixelSpacing,seedPoint,deviation);
        
        
        %% Step 11 Detect Asc Aorta
        
        [ascAorta{iCase}, ascLoc{iCase}, deviation,~] = ascAortaSegmentation(BW,lastDesAorta{iCase},lastDesAortaLoc{iCase},loopIncrementVal,ctpa,pixelSpacing,deviation);
        
        
        %% Step 12 Extract Asc Aorta
        
        [ascDia{iCase},curLoc,preAorta,delBW,seedPoint,deviation,areaList,k,~] = ascAortaExtraction(ascAorta{iCase},ascLoc{iCase},sliceLoc,loopIncrementVal,ctpa,pixelSpacing,deviation);
        
        imgAscaB{iCase} = preAorta ;
        imgAscaO{iCase} = ctpa(:,:,curLoc);
        seedPointOfAscA{iCase} = seedPoint;
        
        ftimeAscA{iCase} = toc(tStartAscA);
                     
        
        %% Step 13 PT Extraction
        
        tStartPT= tic;
        
        [ptLoc{iCase},ptCan{iCase},ptCir{iCase},ptCirD{iCase},isPVDetected,nextLoc,hu{iCase},huN{iCase},preAorta,seedPointPTout,~,~,~,huMask,huMaskLoc] = ptExtraction(PATIENT_POSITION,sliceLoc,curLoc,loopIncrementVal,ctpa,pixelSpacing,seedPoint,deviation,preAorta,delBW,desAorta{iCase},areaList,k);
      
        lastLocPreAorta = ptLoc{iCase};
        seedPointOfAscA{iCase} = getSeedPoint(preAorta);
        
        pt = ptCan{iCase};
        representationOfApicalLevelPVLoc{iCase} = ptLoc{iCase};
        representationOfApicalLevelPV{iCase} = ptCan{iCase};
     
        if(~isPVDetected)
            [ptLoc{iCase},ptCan{iCase},ptCir{iCase},ptCirD{iCase},ptDia{iCase},ptDia3{iCase},ptDiaLines{iCase},ptDiaLines3{iCase},mpa,mpaLoc,~,~,seedPointPTout2, idxOfPre{iCase},indexRho{iCase}] = trackPT(ctpa,ptCan{iCase},nextLoc,loopIncrementVal,pixelSpacing,hu{iCase},preAorta);
            
            seedPointPTout = [seedPointPTout;seedPointPTout2];
            mpaLoc = [nextLoc,mpaLoc];
            mpa = [pt,mpa];
            
            ptSize{iCase} = size(mpaLoc,2);
            
            ptLoc{iCase} = mpaLoc(end);
            temp = mpa(end);
            representationOfApicalLevelPVLoc{iCase} = ptLoc{iCase};
            representationOfApicalLevelPV{iCase} = temp{1};
            
        else
                       
            mpaLoc = ptLoc{iCase};
            ptSize{iCase} = 1;
            [dia,dia3,ptDiaLines{iCase},ptDiaLines3{iCase},idxOfPre{iCase},indexRho{iCase}] =  findOptimalDiameter(pt);
            ptDia{iCase} = round(pixelSpacing*dia);
            ptDia3{iCase} = round(pixelSpacing*dia3);
            mpa = {pt};
            
        end % end of if
        
        ftimePT{iCase} = toc(tStartPT);
        
        
        %% Step 14 HU of DesAorta at Carina Level
        
        if(~isempty(ptLoc{iCase}))
            
            
            currDesALoc = SDdesALoc{iCase};
            currDesA = desAorta{iCase};
            currDeviation = deviation;
            
            desAHU = bsxfun(@times, ctpa(:,:,currDesALoc), cast(~BWt, class(ctpa(:,:,currDesALoc))));
            desAHU = mean2(desAHU(desAHU~=0));
            
            
            %% Curr DesAorta at PT Level
            
            if(size(mpaLoc,2) < 3)
                ptLevel = mpaLoc(1);
            elseif(size(mpaLoc,2) < 5)
                ptLevel = mpaLoc(end-2);
            elseif(size(mpaLoc,2) < 8)
                ptLevel = mpaLoc(end-3);
            else
                ptLevel = mpaLoc(end-4);
            end
            
            disp(strcat(string(iCase)," diff ", string(abs(ptLevel-currDesALoc))))
            
            if (abs(ptLevel-currDesALoc)<30)                
                
                ros = ctpa(:,:,ptLevel);
                
                if(abs(angle) > 25)
                    %ros = imrotate(ros,-angle,'nearest','crop');
                    disp(strcat("abs(angle) > 25 ", string(iCase)))
                end % end of if
                
                ros = ros>(desAHU-round(currDeviation/2)) & ros<(desAHU+round(currDeviation/2));
                
                rosIntersect = ros & currDesA;
                rosIntersect = imfill(rosIntersect,'holes');
                
                C = bwconncomp(rosIntersect);
                stats = regionprops(C,'Area');
                idx =  find([stats.Area]==max([stats.Area]));
                ros = ismember(labelmatrix(C), idx);
                
                [SDdesA{iCase}, desAorta{iCase}, BWt2, seedPoint2, ~] = getDesAortaSD(ros,0,ctpa,ptLevel,pixelSpacing);
                
                imgDesaC{iCase} =  ~BWt2;
                imgDesaB{iCase} = desAorta{iCase};%
                imgDesaO{iCase} = ctpa(:,:,ptLevel);%
                seedPointOfDesA{iCase} = seedPoint2;
                
            end
            
        end % end of if-else
        
        
        %% Display 
        
        %% Trachea
        
        if(showsFigures)
            
            figure,
            RGB = labeloverlay(mat2gray(imgTraO{iCase}),imgTraB{iCase},'Colormap',colorMap,'Transparency',0.25);
            imshow(RGB)
        end      
               
        %% Carina level              
        
        if(showsFigures)
            rMP = rightMidPoint{iCase};
            lMP = leftMidPoint{iCase};                      
            
            figure,
            label_str = ['Distance: ' num2str(distanceCm{iCase},'%0.2f') 'cm'];
            position = [round(rMP(1)) round(rMP(2))-30 abs(round(rMP(1))-round(lMP(1))) 50];
            
            
            RGB = insertObjectAnnotation(mat2gray(imgCarinaO{iCase}),'rectangle',position,label_str,...
                'TextBoxOpacity',0.5,'FontSize',10);
            imshow(RGB)
        end
        %% DesA
        
        if(showsFigures)
            figure,
            RGB = labeloverlay(mat2gray(ctpa(:,:,ptLevel)),desAorta{iCase},'Colormap',colorMap,'Transparency',0.25);
            imshow(RGB)
        end
        
        if(showsFigures)
            figure,
            RGB = insertMarker(mat2gray(ctpa(:,:,ptLevel)),seedPointOfDesA{iCase},'color',colorMap(1,:),'size',5);
            imshow(RGB);
        end
        
        
        %% AscA
        imgAscaB{iCase} = preAorta ;
        imgAscaO{iCase} = ctpa(:,:,lastLocPreAorta);
        seedPointOfAscA{iCase} = seedPoint;
        
        if(showsFigures)
            figure,
            RGB = labeloverlay(mat2gray(ctpa(:,:,lastLocPreAorta)),preAorta,'Colormap',colorMap,'Transparency',0.25);
            imshow(RGB)
        end        
        
        if(showsFigures)
            figure,
            RGB = insertMarker(mat2gray(ctpa(:,:,lastLocPreAorta)),seedPointOfAscA{iCase},'color',colorMap(1,:),'size',5);
            imshow(RGB);
        end
        
        
        %% PT
        
        iPT = round(median(1:size(seedPointPTout,1)));
        
        imgPTB{iCase} = mpa{iPT};
        imgPTO{iCase} = ctpa(:,:,mpaLoc(iPT));%
        
        imgPTDia{iCase} = iPT;        
        
        huMask = reshape(cell2mat(huMask),[512 512 numel(huMask)]);
        humaskNo = round(median(1:size(huMask,3)));
        imgPTHUMask{iCase} =  huMask(:,:,humaskNo);
        imgPTHU{iCase} = ctpa(:,:,huMaskLoc(humaskNo));
        
        if(showsFigures)
            
            figure,
            RGB = labeloverlay(mat2gray(imgPTHU{iCase}),imgPTHUMask{iCase},'Colormap',colorMap,'Transparency',0.48);
            imshow(RGB)
        end
        
        linesP1 = ptDiaLines{iCase};
      
        if(imgPTDia{iCase} > size(linesP1,1))
            linesP =  linesP1(end,:);
        else
            linesP =  linesP1(imgPTDia{iCase},:);
            
        end
        
        if(sum(idxOfPre{iCase}-indexRho{iCase}))
            disp(strcat("failed",string(iCase)));
        end
        
        
        yVal = max([linesP(2);linesP(4)])+12;
        tmpImg = imgPTB{iCase};
        tmpImg(yVal:end,:) =0;
        imgPTB2{iCase} =  tmpImg;          
               
        if(showsFigures)
            figure,imshow(imgPTO{iCase},[])
            hold on
            xy = [linesP(1) linesP(2); linesP(3) linesP(4)];
            plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');
            
            hold off
        end
        
        if(showsFigures)
            figure,
            RGB = labeloverlay(mat2gray(imgPTO{iCase}),imgPTB{iCase},'Colormap',colorMap,'Transparency',0.25);
            imshow(RGB)
        end
        
        seedPointOfPT{iCase} = seedPointPTout(iPT,:);        
  
                
        % apical level of pt/pulmonary valve
        
        imgPVO{iCase} = ctpa(:,:,representationOfApicalLevelPVLoc{iCase});%
        
    
        seedPointOfPV{iCase} = seedPointPTout(end,:); 
        
        counter = counter + 1 ;
        
        tEndMain = toc(tStartMain);
        ftimeMain{iCase} = tEndMain;
        
        %% Displaying Measurments    
        
        SDdesA1 = SDdesA{iCase};
        SDdesA_txt = strcat("Image Noise (SD of HU in a 1 cm² circular ROI in the Descending Aorta) -> ", num2str(SDdesA1), " HU");
       
        ptHU1 = hu{iCase};        
        ptHU_txt = strcat("IV contrast level in Pulmonary Trunk -> ", num2str(ptHU1), " HU");
        ptLoc_s = ptLoc{iCase};  
        
        traDia1 = round(traDia{iCase});
        traDia_txt = strcat("Diameter of the Trachea -> ", num2str(traDia1), " mm");
        
        ascDia1 = round(ascDia{iCase});
        ascDia_txt = strcat("Diameter of the Ascending Aorta -> ", num2str(ascDia1), " mm");
        
        ptDiaF = choosePTDia(ptDia{iCase}); 
        ptDia_txt = strcat("Diameter of the Pulmonary Trunk -> ", num2str(ptDiaF), " mm");        
       
        path_result = fullfile(current_path, "Results/",strcat(file_list(end),'_',timestamp_str,'.txt'));
        
        fileID = fopen(path_result,'w');        
        fprintf(fileID, 'Measurements: \n%s\n%s\n%s\n%s\n%s', SDdesA_txt, ptHU_txt, traDia_txt,ascDia_txt,ptDia_txt);
        fclose(fileID);            

        showResults;       
       
    catch me     
        
        counter = counter + 1 ;        
        fprintf('Error identifier: %s\n  error message: %s\n', me.identifier, me.message);
        disp("Something wrong during the model runtime. The error is logged'");       
                
        path_error = fullfile(current_path, "Logs/",strcat(file_list(end),'_',timestamp_str,'.txt'));
        
        fileID = fopen(path_error,'w');        
        fprintf(fileID, 'Error identifier: %s\n  error message: %s\n', me.identifier, me.message);
        fclose(fileID);
        
    end % end of try-catch
    
end % end of for

toc
tf=cputime-ts;
txt_out = sprintf('The total elapsed time: %d',tf);
fprintf('%s\n', txt_out);














