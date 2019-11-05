FROM ubuntu:xenial

COPY ./*.py /code/
COPY requirements.txt /code/

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
sudo \
wget \
apt-utils \
curl \
git \
dos2unix \
tree \
zip \
unzip \
make \
cmake \
bzip2 \
build-essential \
libtool \
autotools-dev \
automake \
autoconf \
tzdata \
git \
unzip \
ca-certificates \
python3 \
python3-dev \
python3-pip

#=============================================
#install fsl
#=============================================
RUN echo "America/New_York" | sudo tee /etc/timezone && sudo dpkg-reconfigure --frontend noninteractive tzdata
RUN wget -O- http://neuro.debian.net/lists/trusty.de-md.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
#sudo apt-key adv --recv-keys --keyserver pgp.mit.edu 2649A5A9
RUN apt-key adv --recv-keys --keyserver hkp://ha.pool.sks-keyservers.net 0xA5D32F012649A5A9
RUN sudo apt-get update
RUN sudo apt-get install -y fsl=5.0.9-5~nd14.04+1
#=============================================

#=============================================
#install pipeline's python dependencies
#=============================================
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install numpy==1.16.2 setuptools>=40.8.0
RUN python3 -m pip install -r /code/requirements.txt
#=============================================

#=============================================
#install mcr for matlab qsm code
#=============================================
#install mcr for matlab dipole inversion code
ENV MATLAB_VERSION R2016b
RUN mkdir /opt/mcr_install && \
    mkdir /opt/mcr && \
    wget -P /opt/mcr_install https://ssd.mathworks.com/supportfiles/downloads/${MATLAB_VERSION}/deployment_files/${MATLAB_VERSION}/installers/glnxa64/MCR_${MATLAB_VERSION}_glnxa64_installer.zip && \
    unzip -q /opt/mcr_install/MCR_${MATLAB_VERSION}_glnxa64_installer.zip -d /opt/mcr_install && \
    /opt/mcr_install/install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    rm -rf /opt/mcr_install /tmp/*
#=============================================

#mcr
#ENV MCR_VERSION v91
#ENV LD_LIBRARY_PATH /opt/mcr/${MCR_VERSION}/runtime/glnxa64:/opt/mcr/${MCR_VERSION}/bin/glnxa64:/opt/mcr/${MCR_VERSION}/sys/os/glnxa64:/opt/mcr/${MCR_VERSION}/sys/opengl/lib/glnxa64
ENV MCR_INHIBIT_CTF_LOCK 1

#fsl
ENV FSLDIR=/usr/share/fsl/5.0
ENV POSSUMDIR=$FSLDIR
ENV PATH=/usr/lib/fsl/5.0:$PATH
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV FSLMULTIFILEQUIT=TRUE
ENV FSLTCLSH=/usr/bin/tclsh
ENV FSLWISH=/usr/bin/wish
ENV FSLBROWSER=/etc/alternatives/x-www-browser
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0:${LD_LIBRARY_PATH}

WORKDIR /code
RUN chmod +x run.py
ENTRYPOINT ["/code/run.py"]

