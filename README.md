# Automated-Liver-Tumor-Detection

In this project, a computer-aided diagnosis (CAD) system  is provided for Abdominal CT Liver Images.

This implementation comprises of four main phases: 
* liver segmentation
* lesion candidate segmentation
* feature extraction from each candidate lesion
* liver disease classification

Flowchart The proposed method is as follows

![flowchart](https://user-images.githubusercontent.com/87225043/125976327-1487e4ab-031b-481d-9e05-fda150033585.png)

### Preprocessing
In the preprocessing stage, the histogram thresholding method is used to highlight the liver tissue.

<img width="439" alt="4-1" src="https://user-images.githubusercontent.com/87225043/125978787-b332c273-633d-4773-91fc-80c8d536b9b8.png">

### 1.liver segmentation
To liver segmentation, used the levelset method. This step is done in two parts; The first part uses an initial mask to find the area of the liver and finally morphological operations are used to improve the border of the liver image.

<img width="595" alt="4-3" src="https://user-images.githubusercontent.com/87225043/125980032-e645681d-52b7-407e-8645-33092d4355cb.png">


### 2.lesion candidate segmentation
Fast fuzzy c-means (FFCM) clustering is used for lesion candidates extraction.

<img width="552" alt="4-4" src="https://user-images.githubusercontent.com/87225043/125980330-3d28c9e1-9f82-4321-80d7-f1e01ffe9ba1.png">

### 3.feature extraction
In the third phase after extraction of lesion from the liver, variety of features are extracted from each candidate. These features describe the tissue properties and shape of the lesion and include: area, median, mean, skewness, standad deviation, kurtosis and haralick.

### 4.liver disease classification
Finally, these features are used in a classification stage using a support vector machine.
  

