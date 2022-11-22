
function [imgBifOfTrachea,bifurcationPoint,rightMainBronchus,leftMainBronchus,rightMainTracheaSeedPoints,leftMainTracheaSeedPoints,imgAirways,tresholdVal] = airwaysExtraction(BO,PATIENT_POSITION,seedPointOfTrachea,sliceNoOfTrachea)

%AIRWAYSEXTRACTION airways exctraction from the CT-scan
%
%   Examples:
%       [imgBifOfTrachea,bifurcationPoint,rightMainBronchus,leftMainBronchus,rightMainTracheaSeedPoints,leftMainTracheaSeedPoints,imgAirways,tresholdVal] = AIRWAYSEXTRACTION(PATIENT_POSITION,imgMediastinumArea,sliceNoOfAortaSeparation,seedPointOfLeftPA,thresholdValOfAorticArch)

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})


% connected neighborhood parameter for 2D watershed operation
WATERSHED_CONN_VAL = 8;


sliceNumber = size(BO,3);

%% Determine the start order of scan process // top&bottom of ct-scans

if (PATIENT_POSITION == 1)
    lookUpVAlForbifPoint = round(sliceNumber*0.21);
else % idxOfTop == 2 starts with bottom slice
    lookUpVAlForbifPoint = round(sliceNumber - sliceNumber*0.21);
end



%% Step 0 Initial Setup

B2 = BO;
[~,~,TOTAL_SLICE_NO] = size(B2);
% get first and last slice number
[~,loopIncrementVal,lastSliceNo] = getFirstLastSliceLoc(PATIENT_POSITION,TOTAL_SLICE_NO);



%% Step 2 Apply 3D pre-processing for all slices

imgVOI = ~(B2 > -700);

imgTrachea = regiongrowing(imgVOI(:,:,sliceNoOfTrachea),1,[seedPointOfTrachea(2) seedPointOfTrachea(1)]);

%% Step 3 Get trachea


%% Step 4

% pre-allocation
imgAirways = false(size(B2));
imgAirways(:,:,sliceNoOfTrachea) = imgTrachea;


C = bwconncomp(imgTrachea);
statsTra = regionprops(C,'Area','Centroid');
[~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
thisAreaOfTra  = statsTra(indexTra).Area;




%% Step 4.A.2 Smart-Watershed Method
% apply watershed if needed


% Difference between continue and break statement in loop.
% we use continue because we wont apply 3D RG method at the end of algorithm,
% if 3D RG method will be use  then we need to use break instead of continue.

tControl = 71;

%%
%2018-12-21

isPrinted = true;

%%

for iSlice = (sliceNoOfTrachea+loopIncrementVal):loopIncrementVal:lastSliceNo
    
    
    try
        
        if(tControl == iSlice)
            1;
        end
        
        % Bind this parameter to the age of parient later
        if(thisAreaOfTra < 400 )
                      
            if(isPrinted)                
                isPrinted = false;
            end
            if(iSlice == (sliceNoOfTrachea+loopIncrementVal))
                               
                
                imgThresh = (B2(:,:,sliceNoOfTrachea)>-300);
                
                imgCCL = bwconncomp(imgThresh);
                stats = regionprops(imgCCL,'Area','Image');
                idx =  find([stats.Area]==max([stats.Area]));
                imgLung = ismember(labelmatrix(imgCCL), idx);
                
                outerLung = imfill(imgLung,'holes');
                lungArea = outerLung & ~imgThresh;
                imgTemp = bsxfun(@times, B2(:,:,sliceNoOfTrachea), cast( lungArea,class(B2(:,:,sliceNoOfTrachea))));
                
                aveDenOfVOI = sum(sum(imgTemp))/sum(sum(imgTemp<0));
                
                if(aveDenOfVOI > -625)
                    
                    imgVOI = ~(B2(:,:,iSlice) > -400);
                    tresholdVal = -400;
                else
                    imgVOI = ~(B2(:,:,iSlice) > -575);
                    tresholdVal = -575;
                end
                
            else             
                                
                imgThresh = (B2(:,:,iSlice)>-300);
                
                imgCCL = bwconncomp(imgThresh);
                stats = regionprops(imgCCL,'Area','Image');
                idx =  find([stats.Area]==max([stats.Area]));
                imgLung = ismember(labelmatrix(imgCCL), idx);
                
                outerLung = imfill(imgLung,'holes');
                lungArea = outerLung & ~imgThresh;
                imgTemp = bsxfun(@times, B2(:,:,iSlice), cast( lungArea,class(B2(:,:,iSlice))));
                
                aveDenOfVOI = sum(sum(imgTemp))/sum(sum(imgTemp<0));
                
                if(aveDenOfVOI > -625)
                    
                    imgVOI = ~(B2(:,:,iSlice) > -400);
                    tresholdVal = -400;
                else
                    imgVOI = ~(B2(:,:,iSlice) > -575);
                    tresholdVal = -575;
                end
            end
        else
            imgVOI = ~(B2(:,:,iSlice) > -700);
            tresholdVal = -700;
        end
        
        % getting image
        imgVoiFilled = imgVOI;
        %imgVoiFilled = imfill(imgVoiFilled,'holes');
        %imgVoiFilled = imdilate(imgVoiFilled,sse);%check it
        imgForSeedFinding = imgVoiFilled;
        imgForSeedFinding = imgForSeedFinding  & imgTrachea;
        
        
        % get seed point
        thisSeedPointOfTracheaCan = getSeedPoint(imfill(imgForSeedFinding,'holes'));
        
        if(thisSeedPointOfTracheaCan == 0)
            break;
        end % end of if
        
        if(imgVoiFilled(ceil(thisSeedPointOfTracheaCan(2)),ceil(thisSeedPointOfTracheaCan(1))) == 0)
            % % %                 fprintf('RG for imroiFilled missed\n');
            break;
        end % end of if
        
        % in order to imgVoiFilled = imfill(imgVoiFilled,'holes');
        imgVoiFilled = (imgVoiFilled + imgForSeedFinding)>0;
        % calculating new area
        imgRG = regiongrowing(imgVoiFilled,1,[thisSeedPointOfTracheaCan(2) thisSeedPointOfTracheaCan(1)]);
        C = bwconncomp(imgRG);
        statsTra = regionprops(C,'Area','Centroid');
        [~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
        newAreaOfTra = statsTra(indexTra).Area;
        
        
        if(thisAreaOfTra> 400)
            divideFac = 5;
        else
            divideFac = 4;
        end
        
        
        % check calculated area with previous calculation
        % if area is bigger then the previos one then apply watershed
        if (newAreaOfTra < 1750)
            
            %             img1 = imgRG & imgTrachea;
            %             img2 =  (imgTrachea - img1)>0;
            
            img1 = imgRG&imgForSeedFinding;
            img2 = (imgForSeedFinding - img1)>0;
            
            if(max(max(img2)) > 0 )
                img2N = bwlabel(imgVoiFilled);
                img12 = (img2N == max(img2N(img2==1)));
                
                if(~isempty(img12))
                    %if( max(max(abs(imgControl- img12))) ~= 0 )
                    if(sum(sum(abs(imgRG - img12)))>250)% P1C04 == 260
                        
                        % FIXME when use reduceslicethickness this all function
                        % fails at case 01
                        % at least try to see
                        %if(numel(img12(img12==1)) > 1750) % this mean all right or left lung segmented
                        % put here tracking later
                        %end
                        
                        img2 = img12;
                    end
                end
            end
            
            A = bsxfun(@times, B2(:,:,iSlice), cast(img1, 'like', B2(:,:,iSlice)));
            ave1 = mean2(A(A~=0));
            
            A = bsxfun(@times, B2(:,:,iSlice), cast(img2, 'like', B2(:,:,iSlice)));
            ave2 = mean2(A(A~=0));
            
            
            %C = bwconncomp( img1);
            %statsTra = regionprops(C,'Area');
            %[~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
            
            
            C = bwconncomp( img2);
            statsTra = regionprops(C,'Area');
            [~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
            
            if(isempty(indexTra))
                imgTrachea = imgRG;
                imgAirways(:,:,iSlice) =   imgTrachea ;
                
                thisAreaOfTra = newAreaOfTra;
            else
                AreaOfTra = statsTra(indexTra).Area;
                
                if(abs(ave1-ave2) <400)
                    if(AreaOfTra > (thisAreaOfTra/divideFac) && 1750 > AreaOfTra)
                        
                        img2ForSeedFinding = ismember(labelmatrix(C),indexTra(1));
                        img2ForSeedFinding = img2ForSeedFinding & imgVoiFilled;
                           
                        if( max(max(img2ForSeedFinding)) == 0)
                            imgTrachea = imgRG;
                            imgAirways(:,:,iSlice) =   imgTrachea ;
                            
                            thisAreaOfTra = newAreaOfTra;
                        else
                            if(PATIENT_POSITION == 1)
                                temp = iSlice - lookUpVAlForbifPoint;
                            else
                                temp = lookUpVAlForbifPoint - iSlice;
                            end
                            
                            if(temp > 0)
                                bifurcationPoint = iSlice;
                                img2 = (img2 - img1) >0;
                                % get seed point
                                thisSeedPointOfTracheaCan1 = getSeedPoint(imfill(img1,'holes'));
                                thisSeedPointOfTracheaCan2 = getSeedPoint(imfill(img2,'holes'));
                                
                                if(thisSeedPointOfTracheaCan1(1)>thisSeedPointOfTracheaCan2(1))
                                    leftMainTracheaSeedPoints = thisSeedPointOfTracheaCan1;
                                    rightMainTracheaSeedPoints = thisSeedPointOfTracheaCan2;
                                    leftMainBronchus = img1;
                                    rightMainBronchus = img2;
                                else
                                    rightMainTracheaSeedPoints = thisSeedPointOfTracheaCan1;
                                    leftMainTracheaSeedPoints = thisSeedPointOfTracheaCan2;
                                    
                                    leftMainBronchus = img2;
                                    rightMainBronchus = img1;
                                end % end of if-else
                                
                                break;
                            else
                                imgTrachea = imgRG;
                                imgAirways(:,:,iSlice) =   imgTrachea ;
                                
                                thisAreaOfTra = newAreaOfTra;
                                
                            end% end of if-else
                            
                        end % end of if
                        
                    else
                        imgTrachea = imgRG;
                        imgAirways(:,:,iSlice) =   imgTrachea ;
                        
                        thisAreaOfTra = newAreaOfTra;
                    end
                else
                    imgTrachea = imgRG;
                    imgAirways(:,:,iSlice) =   imgTrachea ;
                    
                    thisAreaOfTra = newAreaOfTra;
                end
                
            end % end of if-else
            
        else
            
            imgControl = imgRG;
            % apply watershed
            %imgRG = imfill(imgRG,'holes');
            D = -bwdist(~imgRG,'cityblock');
            D(~imgRG) = -Inf;
            L = watershed(D,WATERSHED_CONN_VAL);
            imgRG(L == 0) = 0;
            
            % FIXME : check if not work, delete it
            % we apply for eroding
            %imgRG = imfill(imgRG,'holes');
            imgVOIWS = imfill(imgRG,'holes');
            %imgVOIWS =  imgVOIWS | (imgVoiFilled  &  imgTrachea); %this one broke p1c34 but save p4c6
            
            [rows, columns] = size(imgTrachea); % Save original size.
            measurements1 = regionprops(imgTrachea, 'Centroid');
            
            imgTracheaSmall = padarray(imgTrachea,[100 100]);
            imgTracheaSmall = imresize(imgTracheaSmall, [rows, columns] ); % Same size as original.
            measurements2 = regionprops(imgTracheaSmall, 'Centroid');
            
            rowsToShift = round(measurements1.Centroid(2)- measurements2.Centroid(2));
            columnsToShift = round(measurements1.Centroid(1) - measurements2.Centroid(1));
            shiftedImage = circshift(imgTracheaSmall, [rowsToShift columnsToShift]);
            
            %this protect current field of trachea
            imgVOIWS =  imgVOIWS | (shiftedImage & imgVoiFilled);
            %imgVOIWS =  imgVOIWS | shiftedImage;
            %imgVOIWS =  imgVOIWS & imgVoiFilled;
            
            % getting image
            imgForSeedFinding = imgVOIWS  &  imgTrachea;
            % regain deleted regions by watershed % seem unnecessary
            %imgForSeedFinding = imgForSeedFinding | imgTrachea;
            
            % getting seed points
            thisSeedPointOfTracheaCan = getSeedPoint(imfill(imgForSeedFinding,'holes'));
            
            if(thisSeedPointOfTracheaCan == 0)
                break;  % FIXME : use "continue", like in segment SC
            end % end of if
            
            if(imgVOIWS(ceil(thisSeedPointOfTracheaCan(2)),ceil(thisSeedPointOfTracheaCan(1))) == 0)
                % % %                     fprintf('RG for imgRG missed\n');
                break;
            end % end of if
            
            % Check area, if bigger tissues were found,re-apply watershed
            % with different parameter
            
            % calculating new area
            % FIXME: SMT looks wrong RG is suitable for this oper.
            
            % combine RG method with logical ope. in order to increase
            % tracheal area
            imgRG = regiongrowing(imgVOIWS,1,[thisSeedPointOfTracheaCan(2) thisSeedPointOfTracheaCan(1)]);
            %imgRG =  imgRG | (imgVoiFilled  &  imgTrachea);
            C = bwconncomp(imgRG);
            statsTra = regionprops(C,'Area');
            [~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
            newAreaOfTra = statsTra(indexTra).Area;
            
            
            if (newAreaOfTra < 1750)
                
                %             img1 = imgRG & imgTrachea;
                %             img2 =  (imgTrachea - img1)>0;
                
                img1 = imgRG&imgForSeedFinding;
                img2 = (imgForSeedFinding - img1)>0;
                
                if(max(max(img2)) > 0 )
                    img2N = bwlabel(imgVoiFilled);
                    img12 = (img2N == max(img2N(img2==1)));
                    
                    if(~isempty(img12))
                        %if( max(max(abs(imgControl- img12))) ~= 0 )
                        if(sum(sum(abs(imgControl - img12)))>250)
                            img2 = img12;
                        end
                    end
                end
                
                A = bsxfun(@times, B2(:,:,iSlice), cast(img1, 'like', B2(:,:,iSlice)));
                ave1 = mean2(A(A~=0));
                
                A = bsxfun(@times, B2(:,:,iSlice), cast(img2, 'like', B2(:,:,iSlice)));
                ave2 = mean2(A(A~=0));
                
                
                C = bwconncomp( img2);
                statsTra = regionprops(C,'Area');
                [~,indexTra] = find([statsTra.Area] == max([statsTra.Area]));
                
                if(isempty(indexTra))
                    imgTrachea = imgRG;
                    imgAirways(:,:,iSlice) =   imgTrachea ;
                    
                    thisAreaOfTra = newAreaOfTra;
                else
                    
                    AreaOfTra = statsTra(indexTra).Area;
                    %regain area
                    
                    if(abs(ave1-ave2) < 400)
                        if(AreaOfTra > (thisAreaOfTra/divideFac) && 1750 > AreaOfTra)
                            
                            img2ForSeedFinding = ismember(labelmatrix(C),indexTra(1));
                            img2ForSeedFinding = img2ForSeedFinding & imgVoiFilled;                            
                                                        
                            if( max(max(img2ForSeedFinding)) == 0)
                                imgTrachea = imgRG;
                                imgAirways(:,:,iSlice) =   imgTrachea ;
                                
                                thisAreaOfTra = newAreaOfTra;
                            else
                                if(PATIENT_POSITION == 1)
                                    temp = iSlice - lookUpVAlForbifPoint;
                                else
                                    temp = lookUpVAlForbifPoint - iSlice;
                                end
                                
                                if(temp > 0)
                                    bifurcationPoint = iSlice;
                                    thisSeedPointOfTracheaCan1 = getSeedPoint(imfill(img1,'holes'));
                                    thisSeedPointOfTracheaCan2 = getSeedPoint(imfill(img2,'holes'));
                                    
                                    if(thisSeedPointOfTracheaCan1(1)>thisSeedPointOfTracheaCan2(1))
                                        leftMainTracheaSeedPoints = thisSeedPointOfTracheaCan1;
                                        rightMainTracheaSeedPoints = thisSeedPointOfTracheaCan2;
                                        
                                        leftMainBronchus = img1;
                                        rightMainBronchus = img2;
                                        
                                    else
                                        rightMainTracheaSeedPoints = thisSeedPointOfTracheaCan1;
                                        leftMainTracheaSeedPoints = thisSeedPointOfTracheaCan2;
                                        
                                        leftMainBronchus = img2;
                                        rightMainBronchus = img1;
                                        
                                    end % end of if-else
                                    
                                    break;
                                else
                                    imgTrachea = imgRG;
                                    imgAirways(:,:,iSlice) =   imgTrachea ;
                                    
                                    thisAreaOfTra = newAreaOfTra;
                                    
                                end% end of if-else
                                
                            end % end of if
                        else
                            imgTrachea = imgRG;
                            imgAirways(:,:,iSlice) =   imgTrachea ;
                            
                            thisAreaOfTra = newAreaOfTra;
                        end
                    else
                        imgTrachea = imgRG;
                        imgAirways(:,:,iSlice) =   imgTrachea ;
                        
                        thisAreaOfTra = newAreaOfTra;
                    end
                    
                end % end of if-else
                
                
            else
                fprintf('Oops! smt goes wrong, watershed does not work!');
            end
            
            
        end % end of if-else
        
    catch Me
        fprintf('An Error occured in trachea segmentation\n');
        disp(Me);
        
    end  % end of try-catch
    
    
end % end of for loop




imgBifOfTrachea = B2(:,:,(bifurcationPoint));

end % end of function


