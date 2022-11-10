# Automated Detection, Segmentation and Measurement of Mediastinal Structures in CT Pulmonary Angiography
#### A fully automated deterministic approach for mediastinal structures analysis 
> Git repo for the manuscript of "Automated Detection, Segmentation and Measurement of Mediastinal Structures in CT Pulmonary Angiography".

> For more detail, please see our [Patent Cooperation Treaty (pct)](https://patents.google.com/patent/WO2022164374A1) application.


# Table of Contents
- [Abstract](#abstract)
- [CADe model](#cade-model)
- [Usage](#usage)


## Abstract

Mediastinal structure measurements are important for the radiologist’s review of computed tomography pulmonary angiography (CTPA) examinations. In the reporting process, radiologists make measurements of diameters, volumes, and organ densities for image quality assessment and risk stratification. However, manual measurement of these features is time consuming. Here, we sought to develop a time-saving automated algorithm that can accurately detect, segment and measure mediastinal structures in routine clinical CTPA examinations. In this study, **700** CTPA examinations collected and annotated. Of these, a training set of **180** examinations were used to develop a fully automated deterministic algorithm. On the test set of **520** examinations, two radiologists validated the detection and segmentation performance quantitatively, and ground truth was annotated to validate the measurement performance. External validation was performed in **47** CTPAs from two independent datasets. The system had **86%-100%** detection and segmentation accuracy in the different tasks. The automatic measurements correlated well to those of the radiologist **(Pearson's r 0.68-0.99)**. Taken together, the fully automated algorithm accurately detected, segmented, and measured mediastinal structures in routine CTPA examinations having an adequate representation of common artifacts and medical conditions.

## CADe Model


![CADe Model](cad_model.png)

**Figure 1. Flowchart of the CADe algorithm.** **(A)** A total of 700 2 mm axial CTPA image stacks were exported from the PACS server. **(B)** In the pre-processing, the linear scale value and the curvature of the CT image were calculated. **(C)** The segmentation chain consisted of four steps, starting with trachea detection followed by DAo, AAo, and PT detection. **(D)** Noise assessment and measurements of mediastinal vascular structures were reported by the system.

## Usage
