# NCBI download environment

```bash
foo@bar$ conda create -n ncbi_download -c conda-forge -c bioconda entrez-direct sra-tools
```
# ABC model


# Juicer

## Get the last stable release
```bash
cd root/5_link_da_enhancers_to_de_genes/src/juicer_1_6
wget https://github.com/aidenlab/juicer/archive/refs/tags/1.6.tar.gz
tar -xf ./1.6.tar.gz
rm ./1.6.tar.gz
mv ./juicer-1.6/* ./juicer
rm -r juicer-1.6/
```

## Get juicer tools jar file and setup juicer
```bash
cd root/5_link_da_enhancers_to_de_genes/src/juicer_1_6/juicer/CPU/common
wget http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar
ln -s juicer_tools.1.7.5_linux_x64_jcuda.0.8.jar juicer_tools.jar
cd root/5_link_da_enhancers_to_de_genes/src/juicer_1_6/
ln -s juicer/CPU scripts
``` 

## Export correct java options 
```bash
export _JAVA_OPTIONS="-Djava.util.prefs.userRoot=/data5/deepro/tmp -Djava.util.prefs.systemRoot=/data5/deepro/tmp -Duser.home=/data5/deepro/tmp -Djava.library.path=/data5/deepro/starrseq/papers/results/5_link_da_enhancers_to_de_genes/src/juicer_1_6/juicer/CPU/common -Xmx128g -Xms49152m"
# remove export java options from juicer.sh script
```

## Create environment to run juicer
```bash
conda create -n hic_pre -c conda-forge -c bioconda openjdk coreutils bwa
```

## Get PPI data from BioGRID
Downloaded from https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.222/BIOGRID-ALL-4.4.222.tab3.zip

