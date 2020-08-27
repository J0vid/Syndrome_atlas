# get shiny serves plus tidyverse packages image
FROM rocker/shiny:3.6.3

# system libraries of general use
RUN apt-get update && apt-get install -y \
    libssl-dev \
    build-essential \
    gfortran \
    libjpeg-dev \
    xorg \
    libx11-dev \
    libglu1-mesa-dev \
    libfreetype6-dev \
    libcurl4-openssl-dev \
    libxml2-dev


# install R packages required
# (change it dependeing on the packages you need)
RUN R -e "install.packages('Morpho', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('geomorph', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('J0vid/Jovid')"
RUN R -e "install.packages('shinycssloaders', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Rvcg', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('visNetwork', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggrepel', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('plotly', repos='http://cran.rstudio.com/')"


# copy the app to the image
COPY data_segmented.Rdata /srv/shiny-server/
COPY app.R /srv/shiny-server/
COPY about_test.html /srv/shiny-server/
COPY test.css /srv/shiny-server/

# select port
EXPOSE 3838

# allow permission
RUN sudo chown -R shiny:shiny /srv/shiny-server

# run app
CMD ["/usr/bin/shiny-server.sh"]
