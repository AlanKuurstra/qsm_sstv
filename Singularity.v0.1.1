Bootstrap: docker
From: ubuntu:xenial

# anaconda2-5.0.1(python 2.7.12)
# nipype 0.13.1
# fsl 5.0.9

#########
%setup
#########
#copy SS_TV matlab scripts
mkdir -p ${SINGULARITY_ROOTFS}/code/matlab_scripts

cp matlab_scripts/SS_TV_script ${SINGULARITY_ROOTFS}/code/matlab_scripts/SS_TV_script
cp matlab_scripts/run_SS_TV_script.sh ${SINGULARITY_ROOTFS}/code/matlab_scripts/run_SS_TV_script.sh

#copy nipype qsm pipeline
cp *.py ${SINGULARITY_ROOTFS}/code/

#########
%post
#########
export DEBIAN_FRONTEND=noninteractive

apt-get update
apt-get install -y --no-install-recommends \
apt-utils \
sudo \
wget \
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
python-fftw \
ca-certificates \
python3-setuptools \
python3-pip

#=============================================
#python3 packages to run r2star script
#=============================================
pip3 install pip==9.0.1
pip3 install numpy==1.13.3 scipy==1.0.0 nibabel==2.2.1
#=============================================

#=============================================
#install fsl
#=============================================
echo "America/New_York" | sudo tee /etc/timezone && sudo dpkg-reconfigure --frontend noninteractive tzdata
wget -O- http://neuro.debian.net/lists/trusty.de-md.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
#sudo apt-key adv --recv-keys --keyserver pgp.mit.edu 2649A5A9
apt-key adv --recv-keys --keyserver hkp://ha.pool.sks-keyservers.net 0xA5D32F012649A5A9
sudo apt-get update
sudo apt-get install -y fsl=5.0.9-4~nd14.04+1
#=============================================

#=============================================
#install anaconda python and nipype
#=============================================
DEST=/opt
mkdir -p $DEST

ANACONDA2_DIR=$DEST/anaconda2
if [ -d $ANACONDA2_DIR ]; then
	rm -rf $ANACONDA2_DIR
fi

#Anaconda2-5.0.1-Linux-x86_64.sh was not compatible with nipype
INST_FILE=Anaconda2-4.2.0-Linux-x86_64.sh 

wget -P $ANACONDA2_DIR --tries=10 https://repo.continuum.io/archive/$INST_FILE
bash $ANACONDA2_DIR/$INST_FILE -b -f -p $ANACONDA2_DIR
rm $ANACONDA2_DIR/$INST_FILE

export PATH=${ANACONDA2_DIR}/bin/:$PATH
conda install -y -c conda-forge nipype=0.13.1
#=============================================

pip install --upgrade pip

#=============================================
#install pygrabbit and pybids
#=============================================
pip install grabbit==0.0.8
#install pybids
git clone https://github.com/INCF/pybids.git /code/pybids
cd /code/pybids
git checkout 37f0bf3f7ce3a1bf6fb31376b8249c324f572048
ln -s /code/pybids/bids /code/bids
#=============================================

#=============================================
#install pyQSM for phase unwrapping and frequency estimation
#=============================================
#install pyQSM for phase unwrapping and frequency estimation
git clone https://github.com/AlanKuurstra/pyQSM.git /code/pyQSM
cd /code/pyQSM
git checkout dfabc133085a4c06f9e4e5643686b7d711da576f
cd /code/pyQSM/unwrap3dInterface
python setup.py build
python setup.py install --install-platlib ..
cd /code/pyQSM/calculateReliabilityInterface
python setup.py build
python setup.py install --install-platlib ..
#=============================================

#=============================================
#install mcr for matlab dipole inversion code
#=============================================
export MATLAB_VERSION=R2016b
mkdir /opt/mcr_install && \
    mkdir /opt/mcr && \
    wget -P /opt/mcr_install http://www.mathworks.com/supportfiles/downloads/${MATLAB_VERSION}/deployment_files/${MATLAB_VERSION}/installers/glnxa64/MCR_${MATLAB_VERSION}_glnxa64_installer.zip && \
    unzip -q /opt/mcr_install/MCR_${MATLAB_VERSION}_glnxa64_installer.zip -d /opt/mcr_install && \
    /opt/mcr_install/install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    rm -rf /opt/mcr_install /tmp/*
#=============================================

#making this dir during setup now
#mkdir -p /code/matlab_scripts 

#make entry point executable
chmod -R +x /code

#########
%files
#########
#initially copied the .py files in this section.
#however, chmod can't run in environment section, it has to run in the post section
#since file copy must happen before post, these copies are now done in setup
#./matlab_scripts/SS_TV_script /code/matlab_scripts/SS_TV_script
#./matlab_scripts/run_SS_TV_script.sh /code/matlab_scripts/run_SS_TV_script.sh
#./*.py /code/

#########
%environment
#########
#tries to chmod when running container, but permission denied
#chmod must be done in post section
#chmod -R +x /code 

#skip all interactive steps
#export DEBIAN_FRONTEND=noninteractive

export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH
export PATH=/opt/anaconda2/bin/:$PATH

#fsl
export FSLDIR=/usr/share/fsl/5.0
export POSSUMDIR=$FSLDIR
export PATH=/usr/lib/fsl/5.0:$PATH
export FSLOUTPUTTYPE=NIFTI_GZ
export FSLMULTIFILEQUIT=TRUE
export FSLTCLSH=/usr/bin/tclsh
export FSLWISH=/usr/bin/wish
export FSLBROWSER=/etc/alternatives/x-www-browser
export LD_LIBRARY_PATH=/usr/lib/fsl/5.0

#mcr
export MCR_VERSION=v91
export LD_LIBRARY_PATH=/opt/mcr/${MCR_VERSION}/runtime/glnxa64:/opt/mcr/${MCR_VERSION}/bin/glnxa64:/opt/mcr/${MCR_VERSION}/sys/os/glnxa64:/opt/mcr/${MCR_VERSION}/sys/opengl/lib/glnxa64:${LD_LIBRARY_PATH}
export=MCR_INHIBIT_CTF_LOCK 1


#########
%runscript
#########
/code/run.py "$@"

