# MorphologyTimeCourses
R Shiny application for time course morphology data

# Timecourse Morpholgy Data Analysis in ***Pseudomonas aeruginosa***  

## Table of contents
1. [Citation](#citation)
2. [Introduction](#content)
3. [Installation](#installation)
4. [Manual](#manual)
5. [Contact](#contact)

----

### Citation <a name="citation"></a>

This R shiny web application has been developed by the Moradigaravand as part of the following working paper:

*Identification of genetic determinants of Pseudomonas aeruginosa biofilm morphology and biofilm formation patterns using a transposon insertion mutant library*

Citation will be added upon the final publication of the manuscript.

----

### Abtsract <a name="content"></a>
Most bacteria are able to adhere to surfaces to grow in complex communities called biofilms. For bacteria growing within biofilms offers protection against environmental stresses, such as antibiotics. This bacterial trait is problematic as it gives human pathogens such as Pseudomonas aeruginosa a way to evade antimicrobial therapy. To understand biofilm formation and to map genes contributing to this process we monitored the biofilm formation using a P. aeruginosa transposon insertion mutant library over time. Hereby we quantified colony size, biofilm color intensity and colony morphology as descriptive biofilm traits. To capture the temporal dynamics, we compared and categorized the patterns of the biofilm traits using a range of timeseries and unsupervised machine learning algorithms. The analysis enabled us to systematically identify genes linked with deviant biofilm formation and morphology patterns. We discovered mutant clusters with distinctive phenotypes of biofilm morphology and colony wrinkleness, encompassing both previously known genes and potential biofilm-forming genes, involved in a wide variety of pathway types. The examination of these genes reveals a group of candidate genes that regulate morphological change in distinct ways over time. We presented our findings as a web server, which represents gene effects on bacterial morphology and can serve as a resource to identify potential treatment targets and/or improve functional annotations. The web application allows the visualization of the growth-related attributes over time, as well as extraacted colonies across four replicates. 


----
### Installation <a name="installation"></a>

There are two ways to run the tool:

- The package may be directly executed as a web application on shinyapps cloud, using the following link:

```
https://daneshmoradigaravand.shinyapps.io/morphology_app/?_ga=2.211242588.854876434.1653884653-409120471.1653027237
```

- The tool is available on DockerHub and may be fetched and run using the following commmands:

```
docker pull daneshmoradigaravand/morphologytimecourses:latest
docker run --rm -p 3838:3838 morphologytimecourses:latest
```

The application is accessible on the following link on the browser:

```
http://0.0.0.0:3838
```

----
### Manual <a name="manual"></a>

The web application has the following feaatures:

<p align="center">
<img src="https://user-images.githubusercontent.com/35295619/170918256-0b3d1073-6455-4405-bd35-6ffcd6456799.jpg" width="500" />
</p>

1- The user chooses the locus tag for the gene of ineterst. 

2- The attributes of the locus tag, including the function, the insertaion site and annotation of the mutant.

3- The comparative table showing the relative value of the quantified features in comparison to other genes. The abnormality score is the relative anomaly score computed from isolatioon forest. 

4- The values of morphology score, biofilm colony size and colony colour intensity over time. The black curve denotes the fitted sigmoid curve from which features in section 3 were extracted.

5- The extracte images of the clones over time. 

----
### Contact <a name="contact"></a>
For queries, please contact [Danesh Moradigaravand](mailto:d.moradigaravand@bham.ac.uk?subject=[GitHub]), Data-Driven Microbiology lab, Center for Computational Biology, University of Birmingham. 
 


-----

