# Dockerized shiny blast server for the pantranscriptome 
# of Astyanax lacustris, a model fish species

# Get rocker/shiny
FROM rocker/shiny

# Create app dir
RUN mkdir /home/shiny-app

# Install the packages
RUN R -e "install.packages(c('XML','plyr','dplyr', 'shinythemes','DT','sfsmisc','rclipboard','shinyjs'))"

# Install Biostrings from Bioconductor
RUN R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install('Biostrings')"

# Copy app to app dir
COPY app.R /home/shiny-app/app.R

# Copy necessary files to app dir
COPY www/ /home/shiny-app/www

# Port to be exposed
# must match the port defined in the app
EXPOSE 8080

# Execute the app
CMD Rscript /home/shiny-app/app.R
