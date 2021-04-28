FROM rocker/r-ver:3.6.5

# ARG WHEN

WORKDIR /IMPC_SP

RUN apt-get update

# Magick++
###RUN apt-get update && apt-get install -y build-essential
###RUN  apt-get install imagemagick
###ADD https://www.imagemagick.org/download/ImageMagick.tar.gz ImageMagick.tar.gz
####RUN mkdir -p ImageMagick
###RUN tar xvzf ImageMagick.tar.gz 
###WORKDIR ImageMagick-7.0.10-27
####RUN cd ImageMagick-7.0.10-27
###RUN echo $(ls -ah)
###RUN chmod 775 configure
###RUN ./configure
###RUN make
###RUN make install 

# XML package
#RUN apt-get -y purge -f libxml2-dev
#RUN apt-get -y clean
#RUN apt-get -y install libxml2 libxml2-dev

# GSL package
#RUN apt install -y libgsl-dev

# zlib
#RUN apt-get -y install zlib1g-dev

# All above
RUN apt-get update && apt-get -y install libxml2 libxml2-dev libgsl-dev zlib1g-dev imagemagick libpq-dev

# data.table R package
# RUN R -e "install.packages('data.table',repos='https://cloud.r-project.org')"
 
# prepare the environment
RUN R -e "source('https://raw.githubusercontent.com/mpi2/impc_stats_pipeline/master/impc_statistical_pipeline/install_dependencies.R')"

# Make the home directory
WORKDIR /IMPC_SP/StatsPipelineReady
# docker run -v ~/mydocker/results:/home/IMPCStatsPipeline  analysis
