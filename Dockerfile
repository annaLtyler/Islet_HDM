FROM rocker/r-base:latest
LABEL maintainer="atyler"
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y libglpk-dev \
	libgmp-dev \
	libxml2-dev \ 
	libcurl4-openssl-dev \
	pandoc \
	pandoc-citeproc

R --slave -e 'install.packages( 'here' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'pheatmap' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'igraph' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'gprofiler2' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'qtl2' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'stringr' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'wordcloud' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'wordcloud2' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'hexbin' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'regress' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'propagate' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'knitr' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'DT' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'RColorBrewer' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'vioplot' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'Matrix' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'RGCCA' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'bnstruct' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'corpcor' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'cluster' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'grid' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'easyPubMed' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'abind' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'GEOquery' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'limma' ,repos='https://cran.rstudio.com/')
R --slave -e 'install.packages( 'pdftools' ,repos='https://cran.rstudio.com/')

R -slave -e 'BiocManager:install("Biobase")'
R -slave -e 'BiocManager:install("BiocGenerics")'
R -slave -e 'BiocManager:install("GEOquery")'
R -slave -e 'BiocManager:install("limma")'

WORKDIR /payload/
CMD ["R"]

