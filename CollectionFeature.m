function [Features, B] = CollectionFeature(Path)
    Files = dir(Path)
    B =[];
    Features = [];
    h = waitbar(0,'Please wait...');
    for j=3:length(Files)
    ParrentPath=sprintf('%s\\%s\\*.jpg',Path,Files(j).name);
    PFiles=dir(ParrentPath);
    for i=1:length(PFiles)
        fn = [ParrentPath(1:end-5) PFiles(i).name];
        LiverImg = {PFiles(i).name}
        categories = {Files(j).name};
        %%%%%%%%%%%%%%%%preprocess%%%%%%%%%%%
        I = imread(fn);  
        I = imresize(I,[256,256]);
        [row col dim]=size(I);
        if dim > 1
           I=rgb2gray(I);
        end

        Pre=imadjust(I,[.4 .7],[0 1]);
        %%%%%%%%%%%%%%%%%%%%%Liver Segmentation%%%%%%%%%%%%%%%%%%%
        F=imadjust(I,[.8 .9],[0 1]);
        z=im2bw(F);
        masklevel = CreateMask(Pre, z);
        I1=imresize(I,[600 600]);
        mask = false(size(I1)); 
        mask(1,1) = true;
        % Compute the weight array based on grayscale intensity differences.
        W = graydiffweight(I1, mask, 'GrayDifferenceCutoff', 25);
        % Segment the image using the weights.
        thresh = 0.01;
        [BW, D] = imsegfmm(W, mask, thresh);
        dd=D(:,:,1)>0.1;
        st=strel('disk',18);

        d1=imerode(dd,st);
        mul=immultiply(d1,I1(:,:,1));
        Img1 = imresize(mul,[256 256]);
        Img=double(Img1(:,:,1));   
        G=fspecial('gaussian',5);
        Img_smooth=conv2(Img,G,'same');  
        [Ix,Iy]=gradient(Img_smooth);
        f1=Ix.^2+Iy.^2;
        g=1./(1+f1);    
        equldis=2; weight=6;  
        mask1=masklevel;
        %[nrow, ncol]=size(I);
        c0=4; 
        initialLSF= -c0*2*(0.5-mask1); 
        u=initialLSF;
        evolution=230;
        % move evolution
        for n=1:evolution
            u=levelset(u, g ,equldis, weight);    
            if mod(n,20)==0
                 pause(1);
                %imshow(Pre, [0, 255]);
                %colormap(gray);hold on;
                %[c,h] = contour(u,[0 0],'r');        
                %hold off;
            end
        end
        
        u=imfill(u,'holes');
        figure;imshow(u);title('u')
        u1=double(imclearborder(im2bw(u)));
        figure;imshow(u1);title('u1')
        g1 = imdilate(u1, strel('disk',3));
        figure;imshow(g1);title('g1')
        g2 = bwmorph(g1, 'open'); 
        figure;imshow(g2);title('g2')
        g2 = uint8(g2);
        liver_img = I.*g2;
        figure;imshow(liver_img);title('liver_img')

        %%%%%%%%%%%%%%%%%%%%%%%%Lesion Segmentation%%%%%%%%%%%%%%%%%%%
        [C,U,LUT,H]=FastFCMeans(liver_img,3); % perform segmentation
        Umap=FM2map(liver_img,U,H);
        BW1 = uint8(Umap(:,:,2));
        BW2 = bwmorph(BW1, 'close'); 
        BW3 = imopen(BW2,strel('disk',3));
        Area = bwarea(BW3);
        BW3=uint8(BW3);

        lesion_img = liver_img.*BW3;
        %%%%%%%%%%%%%%%%%%%%%%%%Feature Extraction%%%%%%%%%%%%%%%%%%%

        Area = bwarea(lesion_img);
        GLCM = graycomatrix(lesion_img,'Offset',[0 1;-1 1;-1 0;-1 -1]);
        stats = grayprops(GLCM, 0);

        Contrast = stats.contr(1);
        Correlation = stats.corrp(1);
        Energy = stats.energ(1);
        Homogeneity = stats.homom1(1);

        Mean = mean(lesion_img(:));
        Standard_Deviation = std2(lesion_img);
        Kurtosis = kurtosis(double(lesion_img(:)));
        Skewness = skewness(double(lesion_img(:)));

        feat = [Contrast,Correlation,Energy,Homogeneity, Area, Mean, Skewness, Standard_Deviation, Kurtosis]



        Features = [Features; feat]
        B=[B;{Files(j).name}]
    end
    end
    close(h)
    Truetype{1,1} = 'Benign';
    Truetype{2,1} = 'Malignant';
    save Truetype Truetype
    %H=[H;{Files(3:end).name}]

end