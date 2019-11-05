# Run QSM Pipeline on Graham

1. Get a compute canada account and ensure you are part of a group with resource allocations

2. ssh into graham:  
ssh **user**@graham.sharcnet.ca

3. Optional: Create a test directory for the example
    <pre>
    cd /scratch/<b>user</b> (eg. /scratch/akuurstr)
    mkdir test
    cd test
    mkdir tar bids qsm
    </pre>

4. Convert data into a tarball:  
    * Get dicoms from cfmm dicom server and convert to tar:  
        singularity run -B **tarLocation**:/output -b **scratchDir**:/scratch **singularityImageLocation**/khanlab-cfmm2tar-master-latest.simg -d **dateSearch** -p **StudySearch** -n **SubjectSearch** /output  
        eg.

        ```
        singularity pull shub://khanlab/cfmm2tar
        singularity run -B tar:/output -B $SCRATCH:/scratch  khanlab-cfmm2tar-master-latest.simg -d 20170112 -p Menon\^Rugby_team -n 2017_01_12_WRT_TBIS1_068 /output
        ```
    * (Alternative) Convert an existing dicom directory to tar:  
        singularity run -B **dicomLocation**:/input -B **tarLocation**:/output **singularityImageLocation**/khanlab-dicom2tar-master-v0.0.3.simg /input /output  
        eg. If you have your dicoms in /scratch/**user**/test/dicom
        
        ```
        singularity pull shub://khanlab/dicom2tar
        singularity run -B dicom:/input -B tar:/output khanlab-dicom2tar-master-latest.simg /input /output
        ```
    

5. Convert tar to bids:  
singularity run -B **tarLocation**:/input -B **bidsLocation**:/output **singularityImageLocation**/khanlab-tar2bids-master-latest.simg -o /output /input/**tarFilename.tar**  
eg. 
```
singularity pull shub://khanlab/tar2bids:0.0.2d
singularity run -B tar:/input -B bids:/output khanlab-tar2bids-master-0.0.2d.simg -o /output /input/Menon_Rugby_team_20170112_2017_01_12_WRT_TBIS1_068_1.9F8E382F.tar
```

6. Download QSM processing bids app:
```
singularity pull shub://AlanKuurstra/qsm_sstv:v0.2.0
```

7. Create a QSM batch program that can be submitted as a job to graham: 
    ```
    vi QSMBatch.sh
    ```
    copy and modify the file:
    <pre>
    #!/bin/bash
    #
    #SBATCH --account=<b>PUT YOUR ACCOUNT HERE (eg. def-akhanf-ab) </b>
    #SBATCH --ntasks 1            # number of tasks
    #SBATCH --cpus-per-task=20
    #SBATCH --mem 64000            # memory pool per process
    #SBATCH -t 1:00:00            # time (D-HH:MM)
    
    export bids_input=<b>full path to bidsLocation (eg. /scratch/akuurstr/test/bids)</b>
    export output=<b>full path to where you want to store your qsm images (eg. /scratch/akuurstr/test/qsm)</b>
    export SINGULARITY_DIR=<b>full path to qsm singularity image (eg. /scratch/akuurstr/test/AlanKuurstra-qsm_sstv-master-v0.2.0.simg)</b>
    
    export subjects=<b>subject list to process (eg. 068)</b>
    
    singularity run \
    -B ${bids_input}:/bids_input \
    -B ${output}:/output \
    ${SINGULARITY_DIR} \
    /bids_input \
    /output \
    participant \
    --participant_label $subjects \
    --SS_TV_lagrange_parameter 0.4 \
    --keep_unnecessary_outputs
    </pre>

8. Submit job:
```
sbatch QSMBatch.sh
```

9. Check queue:
```
squeue -u user 
```




# Run QSM Pipeline Locally
1. Install a virtual machine with ubuntu and singularity by following the steps in the project https://git.cfmm.robarts.ca/cfmm/singularityvagrant.  
Alternative: If your host OS is linux, you can install singularity without a virtual machine:  
   eg. ubuntu:
    <pre>
    sudo wget -O- http://neuro.debian.net/lists/xenial.us-ca.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
    sudo apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9
    sudo apt-get update
    sudo apt-get install singularity-container
    </pre>
2. Follow the above steps 2 - 6 on your local computer.  
**Note: If using cfmm2tar in step 3, you must explicitly indicate a scratch folder instead of using the environment variable $SCRATCH**
3. Instead of steps 7 - 9, run:  
singularity run -B **bidsLocation**:/input -B **qsmLocation**:/output AlanKuurstra-qsm_sstv-master-v0.2.0.simg /input /output participant --participant_label **subjectID** --SS_TV_lagrange_parameter 0.4 --keep_unnecessary_outputs  
eg.

    ```
    singularity run -B bids:/input -B qsm:/output AlanKuurstra-qsm_sstv-master-v0.1.1.simg /input /output participant --participant_label 068 --SS_TV_lagrange_parameter 0.4 --keep_unnecessary_outputs
    ```

### Note that the QSM pipeline can be resource intensive and will fail if it runs out of memory. In this case it is best to use compute canada.
