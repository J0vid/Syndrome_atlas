FROM rocker/shiny
#system libraries of general use
RUN apt-get update && apt-get install -y \
    software-properties-common \
    libssl-dev \
    build-essential


# install R packages required
RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"

#RUN apt-add-repository ppa:marutter/rrutter3.5
#RUN apt-add-repository ppa:marutter/c2d4u3.5
RUN apt-add-repository ppa:zarquon42/statismo-develop
RUN sudo apt update --allow-insecure-repositories
RUN sudo apt install -y \
    statismo-dev    
# install R packages required
RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"
#RUN R -e "devtools::install_github("zarquon42b/RvtkStatismo",ref="develop")"
RUN R -e "devtools::install_github('dmurdoch/rgl')"
RUN R -e "install.packages('future', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('promises', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Morpho', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('geomorph', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('j0vid/Jovid')"
RUN R -e "install.packages('Rvcg', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('sparsediscrim', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyBS', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinycssloaders', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('plotly', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggrepel', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('grid', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('visNetwork', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinydashboard', repos='http://cran.rstudio.com/')"
