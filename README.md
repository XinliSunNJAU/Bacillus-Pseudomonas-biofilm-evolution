Population Sequencing Mutations Table Analysis
Project Description
This R script performs a detailed analysis of genetic mutations in a evolution study. It includes steps for loading necessary packages, reading mutation data from an Excel file, filtering mutations based on specific criteria, and summarizing the mutations. The script ends with creating various visualizations to understand mutation patterns across different lineages and times.

Installation
Prerequisites
R
RStudio (Recommended)
Required Packages
This script requires the following R packages:

tidyverse
readxl
openxlsx
ggplot2
ggpubr
ggsci
vegan
ape
gridExtra
cowplot
multcompView
agricolae
reshape2
lmerTest
brms
emmeans

Usage
Run the script in RStudio or any R environment. Ensure the data file Mutation.all.xlsx is in your working directory, or modify the path in the script to match your data's location.

Features
Automatically installs required packages.
Filters mutation data based on ancestral status, frequency, and type.
Summarizes mutations and visualizes data for different lineages and time points.
Contributing
Contributions to this project are welcome. Please fork the repository, make your changes, and submit a pull request.

License
This project is licensed under the MIT License - see the LICENSE.md file for details.

Contact
For any questions or issues, please open an issue on the GitHub repository page or contact sunxinli@njau.edu.cn.

Acknowledgments
Thanks to all the contributors who have helped with testing, feedback, and suggestions.
