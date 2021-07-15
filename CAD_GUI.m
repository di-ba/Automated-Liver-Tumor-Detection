function varargout = CAD_GUI(varargin)
% CAD_GUI MATLAB code for CAD_GUI.fig
%      CAD_GUI, by itself, creates a new CAD_GUI or raises the existing
%      singleton*.
%
%      H = CAD_GUI returns the handle to a new CAD_GUI or the handle to
%      the existing singleton*.
%
%      CAD_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAD_GUI.M with the given input arguments.
%
%      CAD_GUI('Property','Value',...) creates a new CAD_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CAD_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CAD_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CAD_GUI

% Last Modified by GUIDE v2.5 11-Jun-2019 13:13:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CAD_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CAD_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CAD_GUI is made visible.
function CAD_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CAD_GUI (see VARARGIN)

% Choose default command line output for CAD_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CAD_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CAD_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.jpg;*.png;*.bmp','Pick an CT Image');
if isequal(FileName,0)||isequal(PathName,0)
    warndlg('User Press Cancel');
else
    I = imread([PathName,FileName]);
    I = imresize(I,[256,256]);
    [row col dim]=size(I);
if dim > 1
   I=rgb2gray(I);
end
  
  axes(handles.axes1)
  imshow(I);%title('CT Image');
 
  handles.ImgData = I;
  
  guidata(hObject,handles);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'ImgData')
    I = handles.ImgData;
    
    Pre=imadjust(I,[.4 .7],[0 1]);
    axes(handles.axes2)
    imshow(Pre);%title('CT Image');

    handles.ImgPre = Pre;

guidata(hObject,handles);
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = handles.ImgData;
F=imadjust(I,[.8 .9],[0 1]);
Pre = handles.ImgPre;
z=im2bw(F);

I2=I;

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

% Img1 = imresize(mul,[256 256]);
% Img=double(Img1(:,:,1));   
% G=fspecial('gaussian',5);
% Img_smooth=conv2(Img,G,'same');  
% [Ix,Iy]=gradient(Img_smooth);
% f1=Ix.^2+Iy.^2;
% g=1./(1+f1);    
% equldis=2; weight=6; 
% 
% 
% mask1=masklevel;
% %[nrow, ncol]=size(I);
% c0=4; 
% initialLSF= -c0*2*(0.5-mask1); 
% u=initialLSF;
% evolution=230;


Img1 = imresize(mul,[256 256]);
Img=double(Img1(:,:,1));   
G=fspecial('gaussian',5);
Img_smooth=conv2(Img,G,'same');  
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);    
equldis=2; weight=6;   
width = 256;
height = 256;
radius = 10;
centerW = width/3.3;
centerH = height/2.3;
[W,H] = meshgrid(1:width,1:height);
mask = ((W-centerW).^2 + (H-centerH).^2) < radius^2;


%  mask=roipoly(Img1)
if  mean2(I2)>50
mask=imread('mask1.jpg');
else
mask=imread('mask.jpg');
end

BW = double(im2bw(mask)); 
% BW=mask;
[nrow, ncol]=size(Img1);
c0=4; 
initialLSF= -c0*2*(0.5-BW); 
u=initialLSF;
u=initialLSF;
evolution=230;



% move evolution
for n=1:evolution
    u=levelset(u, g ,equldis, weight);    
    if mod(n,20)==0
         pause(1);
        axes(handles.axes3)
        imshow(I, [0, 255]);colormap(gray);hold on;
        [c,h] = contour(u,[0 0],'r'); 
        w = h.LineWidth;
        h.LineWidth = 2;
        hold off;
    end
end

u=imfill(u,'holes');
u1=double(imclearborder(im2bw(u)));
g2 = imdilate(u1, strel('disk',3));
BW2 = bwmorph(g2, 'open'); 
BW2 = uint8(BW2);
img_out = I.*BW2;

axes(handles.axes3)
imshow(img_out);

handles.ImgLiver = img_out;

BW_groundTruth= imread('Images\Malignant\i_57b.jpg');
% gt=BW_groundTruth(:);
% segm=BW_Created(:);
BW_groundTruth = logical(BW_groundTruth);
BW_Created = logical(BW2);
similarityDICE = dice(BW_Created, BW_groundTruth)
similarityJAC = jaccard(BW_Created, BW_groundTruth)
%similarityFSCORE = bfscore(BW_Created, BW_groundTruth)
corr= corr2(BW_Created, BW_groundTruth)


guidata(hObject,handles);




% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
liver_img  = handles.ImgLiver;
[C,U,LUT,H]=FastFCMeans(liver_img,3); % perform segmentation

Umap=FM2map(liver_img,U,H);

BW1 = uint8(Umap(:,:,2));

BW2 = bwmorph(BW1, 'close'); 
%BW3 = bwmorph(BW2, 'open'); 
BW3 = imopen(BW2,strel('disk',3));
Area = bwarea(BW3);
BW3=uint8(BW3);

lesion_img = liver_img.*BW3;

axes(handles.axes4)
imshow(lesion_img);

handles.ImgLesion = lesion_img;

guidata(hObject,handles);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lesion_img = handles.ImgLesion;
whos lesion_img
Area = bwarea(lesion_img)
GLCM = graycomatrix(lesion_img,'Offset',[0 1;-1 1;-1 0;-1 -1]);
stats = grayprops(GLCM, 0);

Contrast = stats.contr(1)
Correlation = stats.corrp(1);
Energy = stats.energ(1);
Homogeneity = stats.homom1(1);

Mean = mean(lesion_img(:))
Standard_Deviation = std2(lesion_img);
%Standard_Deviation1 = std(lesion_img(:));
Kurtosis = kurtosis(double(lesion_img(:)));
Skewness = skewness(double(lesion_img(:)));

feat = [Contrast,Correlation,Energy,Homogeneity, Area, Mean, Skewness, Standard_Deviation, Kurtosis]

set(handles.uitable1,'Data',feat);
set(handles.uitable1, 'ColumnName', {'Contrast', 'Correlation','Energy','Homogeneity',....
        'Area','Mean','Skewness','Standard_Deviation','Kurtosis'});
set(handles.uitable1, 'RowName', {'Value'});
handles.Feature = feat;
guidata(hObject,handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
testFeature = handles.Feature
load TrainFeatures
load Truetype
%******************Feature Selection*************************
X = TrainImgFeatures;
Y = Livercate';

 svmStruct1 = svmtrain(X,Y,'kernel_function', 'rbf');
 species = svmclassify(svmStruct1,testFeature)
 s = species{1};
     if isequal(s,Truetype{1})
         Imgcate=species{1}
         helpdlg(' Benign Tumor ');
         disp(' Benign Tumor ');
     elseif isequal(s,Truetype{2})

         Imgcate = species{1}
         helpdlg(' Malignant Tumor ');
         disp(' Malignant Tumor ');
     end

 set(handles.edit1,'string',Imgcate);
 guidata(hObject,handles);



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Path = handles.Path;
[TrainImgFeatures, Livercate] = CollectionFeature(Path);

save TrainFeatures 'TrainImgFeatures' 'Livercate'
guidata(hObject,handles);

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

folder_name = uigetdir;
if isequal(folder_name,0)
    warndlg('User Press Cancel');
else
    handles.Path = folder_name;
end
guidata(hObject,handles);
