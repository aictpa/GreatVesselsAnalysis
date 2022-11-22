% SHOWRESULTS  Illustration of measurments
%

%   Copyright 2022
%   Author  - Kahraman A. Teymur
%   Email   - ali_teymur*kahraman=igp*uu*se
%   Real_email = regexprep(Email,{'=','*'},{'@','.'})

txt_version = version('-release');

version_year = str2double(txt_version(1:4));
version_period = lower(txt_version(5));


%% Trachea

label_str = ['Trachea Diameter (transverse): ' num2str(traDia1(iCase),'%0.2f') ' mm'];

props = regionprops(imgTraB{iCase}, 'Area', 'Perimeter','Centroid','EquivDiameter','BoundingBox');
centroids = round(cat(1, props.Centroid));
radius = floor([props.EquivDiameter]/2);
linesA = [(centroids(1)-radius) (centroids(2)) (centroids(1)+radius) (centroids(2))];

if version_year > 2020 || (version_year == 2020 && version_period == 'b' )    
    RGB = insertText(mat2gray(imgTraO{iCase}),[0 0],label_str,'FontSize',14,'BoxColor',...
        'red','BoxOpacity',0.5,'TextColor','white');
    RGB = insertObjectMask(RGB, imgTraB{iCase},'LineColor','red','LineWidth',1,'Opacity',0.05);
    figure,imshow(RGB,[])
    hold on
    xyA = [linesA(1) linesA(2); linesA(3) linesA(4)];
    plot(xyA(:,1),xyA(:,2),'LineWidth',2,'Color','blue');
    hold off
    
else
    position = ceil(props.BoundingBox);
    RGB = insertObjectAnnotation(mat2gray(imgTraO{iCase}),'rectangle', position,label_str,...
        'TextBoxOpacity',0.5,'FontSize',10);
    figure,imshow(RGB)
    hold on
    xyA = [linesA(1) linesA(2); linesA(3) linesA(4)];
    plot(xyA(:,1),xyA(:,2),'LineWidth',2,'Color','red');
    hold off
    
end


%% AAo

label_str = ['AAo Diameter: ' num2str(ascDia1(iCase),'%0.2f') ' mm'];

props = regionprops(imgAscaB{iCase}, 'Area', 'Perimeter','Centroid','EquivDiameter','BoundingBox');
centroids = round(cat(1, props.Centroid));
radius = round([props.EquivDiameter]/2);
linesA = [(centroids(1)-radius) (centroids(2)) (centroids(1)+radius) (centroids(2))];

if version_year > 2020 || (version_year == 2020 && version_period == 'b' )    
    RGB = insertText(mat2gray(imgAscaO{iCase}),[0 0],label_str,'FontSize',14,'BoxColor',...
        'red','BoxOpacity',0.5,'TextColor','white');
    RGB = insertObjectMask(RGB, imgAscaB{iCase},'LineColor','red','LineWidth',1,'Opacity',0.05);
    figure,imshow(RGB,[])
    hold on
    xyA = [linesA(1) linesA(2); linesA(3) linesA(4)];
    plot(xyA(:,1),xyA(:,2),'LineWidth',2,'Color','blue');
    hold off
    
else
    position = round(props.BoundingBox);
    RGB = insertObjectAnnotation(mat2gray(imgAscaO{iCase}),'rectangle', position,label_str,...
        'TextBoxOpacity',0.5,'FontSize',10);
    figure,imshow(RGB)
    hold on
    xyA = [linesA(1) linesA(2); linesA(3) linesA(4)];
    plot(xyA(:,1),xyA(:,2),'LineWidth',2,'Color','red');
    hold off
    
end

    
%% PT

label_str = ['PT Diameter: ' num2str(ptDiaF(iCase),'%0.2f') ' mm'];
    
if version_year > 2020 || (version_year == 2020 && version_period == 'b' )    
    RGB = insertText(mat2gray(imgPTO{iCase}),[0 0],label_str,'FontSize',14,'BoxColor',...
    'red','BoxOpacity',0.5,'TextColor','white');
    RGB = insertObjectMask(RGB,imgPTB{iCase},'LineColor','red','LineWidth',1,'Opacity',0.05);
    figure,imshow(RGB,[])
    hold on
    xy = [linesP(1) linesP(2); linesP(3) linesP(4)];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','blue');
    hold off
else
    stats = regionprops(imgPTB{iCase},'BoundingBox');
    position = round(stats.BoundingBox);
    RGB = insertObjectAnnotation(mat2gray(imgPTO{iCase}),'rectangle', position,label_str,...
        'TextBoxOpacity',0.5,'FontSize',10);
    figure,imshow(RGB)
    hold on
    xy = [linesP(1) linesP(2); linesP(3) linesP(4)];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');
    hold off
    
end


%% IV PT

label_str = ['IV Contrast in PT: ' num2str(ptHU1(iCase),'%0.f') ' HU'];

if version_year > 2020 || (version_year == 2020 && version_period == 'b' )   
    RGB = insertText(mat2gray(imgPTHU{iCase}),[0 0],label_str,'FontSize',14,'BoxColor',...
    'red','BoxOpacity',0.5,'TextColor','white');
    RGB = insertObjectMask(RGB, imgPTHUMask{iCase},'LineColor','red','LineWidth',1,'Opacity',0.05);
    figure,imshow(RGB,[])
else
    stats = regionprops(imgPTHUMask{iCase},'Centroid','MajorAxisLength','MinorAxisLength');    
    centers = floor(stats.Centroid);
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = floor(diameters/2);
    position = [centers radii];
    RGB = insertObjectAnnotation(mat2gray(imgPTHU{iCase}),'circle', position,label_str,...
        'TextBoxOpacity',0.5,'FontSize',10);
    figure,imshow(RGB) 
end


%% Image Noise

label_str = ['Image Noise: ' num2str(SDdesA1(iCase),'%0.f') ' HU'];

if version_year > 2020 || (version_year == 2020 && version_period == 'b' )
    
    RGB = insertText(mat2gray(imgDesaO{iCase}),[0 0],label_str,'FontSize',14,'BoxColor',...
    'red','BoxOpacity',0.5,'TextColor','white');
    RGB = insertObjectMask(RGB,  imgDesaC{iCase},'LineColor','red','LineWidth',1,'Opacity',0.05);
    figure,imshow(RGB,[])
    
else  
    stats = regionprops( imgDesaC{iCase},'Centroid','MajorAxisLength','MinorAxisLength');    
    centers = floor(stats.Centroid);
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    radii = floor(diameters/2);
    position = [centers radii];
    RGB = insertObjectAnnotation(mat2gray(imgDesaO{iCase}),'circle', position,label_str,...
        'TextBoxOpacity',0.5,'FontSize',10);
    figure,imshow(RGB)     
end




