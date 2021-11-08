library(ggplot2)
#MS2 threshold
ms2thr = 0.01
ms2_per = 0.5
ms2_n = 5
#MS1 threshold, peak height
sthr = 10000
ms2count_per = 0.05
#directory
dir_path = "E:/hrms"

# read the ms2 library
ms2_sj = t(read.csv(file.choose(), stringsAsFactors = F, sep = "\t", quote = ""))
# read the feature table
ms_data = read.csv(file.choose(), stringsAsFactors = F, sep = "\t", quote = "")

#==================extract key/assistance ions=================================
ms2_gt = data.frame(matrix(rep(0,30000), nrow = 100, ncol = 300))
rown = 0
coln = 0 
ms2col_l = vector(mode = "numeric")

jrow = vector(mode = "numeric")
jcol = vector(mode = "numeric")

for (wz in 1:nrow(ms2_sj)){
  rown = rown + 1
  for (lz in 1:ncol(ms2_sj)){
    if (!is.na(ms2_sj[wz,lz])){
      
      if (length(ms2col_l[abs(ms2col_l - ms2_sj[wz,lz]) < ms2thr]) == 0){
        coln = length(ms2col_l) + 1
        ms2col_l = c(ms2col_l,ms2_sj[wz,lz])
        ms2_gt[rown, coln] = 1
      }
      else{
        coln = which(abs(ms2col_l - ms2_sj[wz,lz]) < ms2thr)[1]
        ms2_gt[rown, coln] = 1
      }
      
    }
  }
}

for (i in 1:ncol(ms2_gt)){
  if (sum(ms2_gt[,i]) <= 1){
    jcol = c(jcol,i)
  }
}

for (i in 1:nrow(ms2_gt)){
  if (sum(ms2_gt[i,]) <= 1){
    jrow = c(jrow,i)
  }
}

ms2_gt = ms2_gt[-jrow,]
ms2_gt = ms2_gt[,-jcol]
ms2col_l = ms2col_l[-jcol]

ms2_kl = vector(mode = "numeric")
ms2_al = vector(mode = "numeric")
ms2_sk = vector(mode = "numeric")
ms2_sa = vector(mode = "numeric")
ms2_sdf = data.frame(kl = vector(), al = vector())
ms2_ldf = data.frame(kl = vector(), al = vector())

for (i in 1:ncol(ms2_gt)){
  if (sum(ms2_gt[,i]) >= (nrow(ms2_sj) * ms2_per)){
    print(sum(ms2_gt[,i]))
    print(ms2col_l[i])
    ms2_kl = c(ms2_kl, ms2col_l[i])
    ms2_sk = c(ms2_sk,i)
  }
  if (sum(ms2_gt[,i]) > ms2_n && sum(ms2_gt[,i]) < nrow(ms2_sj) * ms2_per){
    print(sum(ms2_gt[,i]))
    print(ms2col_l[i])
    ms2_al = c(ms2_al, ms2col_l[i])
    ms2_sa = c(ms2_sa, i)
  }
}

for (i in 1:nrow(ms2_sj)){
  lklz = vector(mode = "numeric")
  sklz = vector(mode = "numeric")
  lalz = vector(mode = "numeric")
  salz = vector(mode = "numeric")
  for (lz in ms2_sk){
    if (ms2_gt[i,lz] == 1){
      lklz = c(lklz, ms2col_l[lz])
      sklz = c(sklz,lz)
    }
  }
  for (lz in ms2_sa){
    if (ms2_gt[i,lz] == 1){
      lalz = c(lalz, ms2col_l[lz])
      salz = c(salz,lz)
    }
  }
  ms2_sdf = rbind(ms2_sdf, data.frame(kl = I(list(sklz)), al = I(list(salz))))
  ms2_ldf = rbind(ms2_ldf, data.frame(kl = I(list(lklz)), al = I(list(lalz))))
}
rownames(ms2_sdf) = rownames(ms2_sj)
rownames(ms2_ldf) = rownames(ms2_sj)

#======================================================================


ms2c <- function(msstr) {
  
  ms2_list = unlist(strsplit(msstr, " "))
  
  ms2_int_c = vector(mode = "numeric")
  
  for (i in 1:length(ms2_list)){
    ms2mzc = unlist(strsplit(ms2_list[i],":"))[2]
    ms2_int_c = c(ms2_int_c,as.numeric(ms2mzc))
  }
  
  ms2_mc = max(ms2_int_c) * ms2count_per
  return(ms2_mc)
}

ms2zh <- function(msstr, ms2int) {
  
  wzsk = vector(mode = "numeric")
  wzsa = vector(mode = "numeric")
  
  ms2_list = unlist(strsplit(msstr, " "))
  for (lz in ms2_list){
    ms2mz = as.numeric(unlist(strsplit(lz, ":"))[1])
    ms2mzc = as.numeric(unlist(strsplit(lz, ":"))[2])
    
    if (ms2mzc > ms2int && length(ms2_kl[abs(ms2mz - ms2_kl) < ms2thr]) > 0){
      index = which(abs(ms2mz - ms2col_l) < ms2thr)[1]
      wzsk = c(wzsk, index)
    }
    if (ms2mzc > ms2int && length(ms2_al[abs(ms2mz - ms2_al) < ms2thr]) > 0){
      index = which(abs(ms2mz - ms2col_l) < ms2thr)[1]
      wzsa = c(wzsa, index)
    }
    
  }
  return(list(wzsk,wzsa))
}

ms2zhd <- function(msstr) {
  
  wzsk = vector(mode = "numeric")
  wzsa = vector(mode = "numeric")
  
  for (lz in msstr){
    ms2mz = as.numeric(lz)
    
    if (length(ms2_kl[abs(ms2mz - ms2_kl) < ms2thr]) > 0){
      index = which(abs(ms2mz - ms2col_l) < ms2thr)[1]
      wzsk = c(wzsk, index)
    }
    if (length(ms2_al[abs(ms2mz - ms2_al) < ms2thr]) > 0){
      index = which(abs(ms2mz - ms2col_l) < ms2thr)[1]
      wzsa = c(wzsa, index)
    }
    
  }
  return(list(wzsk,wzsa))
}

jac <- function(ms2s,rs) {
  ms2s = unlist(ms2s)
  rs = unlist(rs)
  its = length(intersect(ms2s,rs))
  uts = length(union(ms2s,rs))
  return(its / uts)
}

fs_thr <- function(msd, kjacd, ajacd) {
  kr = vector(mode = "numeric")
  ar = vector(mode = "numeric")
  wza = vector(mode = "numeric")
  for (wz in 1:nrow(msd)){
    if (ms_data[wz,8] > sthr){
      ms2_count = ms2c(ms_data[wz,31])
      sp = ms2zh(msd[wz,31], ms2_count)
      for (nr in 1:nrow(ms2_sdf)){
        kjac = jac(sp[[1]], unclass(ms2_sdf[nr,1]))
        kr = c(kr, kjac)
        ajac = jac(sp[[2]], unclass(ms2_sdf[nr,2]))
        ar = c(ar, ajac)
        if (kjac > kjacd && ajac > ajacd && !(wz %in% wza)) {
          wza = c(wza, wz)
          print(wz)
          print(paste(kjac, ajac, " "))
        }
      }
    }
  }
  return(wza)
}

scr_thr <- function(msarr) {
  
  for (wz in msarr){
    num = nrow(ms2_sdf)
    lzdf = data.frame(kl = numeric(num), al = numeric(num), kah = numeric(num), zl = character(num), stringsAsFactors = F) 
      
      ms2_count = ms2c(ms_data[wz,31])
      sp = ms2zh(ms_data[wz,31], ms2_count)
      
      for (nr in 1:nrow(ms2_sdf)){
        kjac = jac(sp[[1]], unclass(ms2_sdf[nr,1]))
        ajac = jac(sp[[2]], unclass(ms2_sdf[nr,2]))
        lzdf[nr, 1] = ifelse(kjac > 0, kjac, 0)
        lzdf[nr, 2] = ifelse(ajac > 0, ajac, 0)
        lzdf[nr, 3] = 1.5*kjac + ajac 
        lzdf[nr, 4] = rownames(ms2_sdf)[nr]
      }
      
      lzdf = lzdf[order(lzdf[,3]), ]
      lzdfp = data.frame(lzc = c(lzdf[26:30, 1],lzdf[26:30, 2]), lzf = c(rep("kl", 5), rep("al", 5)), zl = lzdf[26:30,4])
      fn = paste(wz,".jpeg",sep = "")
      jpeg(fn)
      p = ggplot(data = lzdfp, mapping = aes(x = factor(zl,levels = as.factor(lzdf[26:30, 4])), y = lzc, fill = lzf))  + 
        geom_bar(stat="identity", position = "dodge") + labs(x = "wz", y = "jac", fill = "lz") + theme_bw() +
        theme(axis.text.x = element_text(angle=25, size = 10, vjust = 0.6))
      print(p)
      dev.off()
    }
}


#fast screening(kjac score1 ajac score2)
res = fs_thr(ms_data, kjac, ajac)
# visualize the results
scr_thr(res)