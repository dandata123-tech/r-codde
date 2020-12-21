options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options("repos" = c(CRAN="http://mirrors.cloud.tencent.com/CRAN/")) 
options(download.file.method = 'libcurl')
options(url.method='libcurl')

BiocManager::install("miRNAtap",ask = F,update = F)
BiocManager::install("topGO",ask = F,update = F)
BiocManager::install("miRNAtap.db",ask = F,update = F)