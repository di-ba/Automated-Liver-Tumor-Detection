function u = LiverSegmentation(I, masklevel)
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
            axes(handles.axes3)
            imshow(Pre, [0, 255]);colormap(gray);hold on;
            [c,h] = contour(u,[0 0],'r');        
            hold off;
        end
    end


end