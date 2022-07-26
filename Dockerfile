FROM rocker/r-ver:4.1.2
RUN apt-get update && apt-get install -y  cmake git-core libcurl4-openssl-dev libgit2-dev libicu-dev libpng-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc zlib1g-dev && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN R -e 'install.packages("remotes")'
RUN Rscript -e 'remotes::install_version("rlang",upgrade="never", version = "1.0.2")'
RUN Rscript -e 'remotes::install_version("pkgload",upgrade="never", version = "1.2.4")'
RUN Rscript -e 'remotes::install_version("testthat",upgrade="never", version = "3.1.2")'
RUN Rscript -e 'remotes::install_version("abind",upgrade="never", version = "1.4-5")'
RUN Rscript -e 'remotes::install_version("ggplot2",upgrade="never", version = "3.3.5")'
RUN Rscript -e 'remotes::install_version("dplyr",upgrade="never", version = "1.0.8")'
RUN Rscript -e 'remotes::install_version("tidyr",upgrade="never", version = "1.2.0")'
RUN Rscript -e 'remotes::install_version("vroom",upgrade="never", version = "1.5.5")'
RUN Rscript -e 'remotes::install_version("shiny",upgrade="never", version = "1.7.1")'
RUN Rscript -e 'remotes::install_version("config",upgrade="never", version = "0.3.1")'
RUN Rscript -e 'remotes::install_version("shinytest")'
RUN Rscript -e 'remotes::install_version("vdiffr",upgrade="never", version = "1.0.4")'
RUN Rscript -e 'remotes::install_version("markdown",upgrade="never", version = "1.1")'
RUN Rscript -e 'remotes::install_version("reshape2",upgrade="never", version = "1.4.4")'
RUN Rscript -e 'remotes::install_version("plotly",upgrade="never", version = "4.10.0")'
RUN Rscript -e 'remotes::install_version("ggpubr",upgrade="never", version = "0.4.0")'
RUN Rscript -e 'remotes::install_version("DT",upgrade="never", version = "0.21")'
RUN Rscript -e 'remotes::install_version("JBrowseR",upgrade="never", version = "0.9.0")'
RUN Rscript -e 'remotes::install_version("golem",upgrade="never", version = "0.3.2")'
RUN Rscript -e 'remotes::install_version("shinydashboard",upgrade="never", version = "0.7.2")'
RUN Rscript -e 'remotes::install_version("shinyWidgets",upgrade="never", version = "0.7.0")'
RUN Rscript -e 'remotes::install_version("shinythemes",upgrade="never", version = "1.2.0")'
RUN Rscript -e 'remotes::install_version("shinyjs",upgrade="never", version = "2.0.0")'
RUN mkdir /build_zone
ADD . /build_zone
WORKDIR /build_zone
RUN R -e 'remotes::install_local(upgrade="never")'
RUN rm -rf /build_zone
EXPOSE 80
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');viewpoly::run_app()"
