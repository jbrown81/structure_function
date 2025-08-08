# Copyright (c) 2025 Radiata, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

library(mgcv)
library(gratia)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)

tests=c("MMSETot", "TRAILB", "MTTime", "StrpCor", "StrpCNCor", "TRAILA", "DFCorr", "DCorr", "BNTTot", "ANIMALS", "VEG", "PPVTtotal", "Verbal", "Syntax", "DigitFW", "DigitBW", "ModRey", "NumbLoc", "Rey10m", "TrCoTot_Short", "Corr10", "Corr30", "RecogFP", "cdr_ftld_box", "MMSETot_raw")
test_flip_sign=c(1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1)

all_models=list()
model_devs <- data.frame()
term_ts <- data.frame()
term_t_ps <- data.frame()
term_fs <- data.frame()
term_f_ps <- data.frame()
term_f_dfs <- data.frame()
all_preds=list()
for (i in 1:25) {
  cur_file=sprintf("/Users/jbrown2/Desktop/adftd_structure_function/cognition_data/yX_%s.txt",tests[i])
  cur_data<-read.table(cur_file)
  colnames(cur_data) <- c("cog","s1","f1","s1xf1","s2","f2","s2xf2","s3","f3","s3xf3","age","sex","scanner","motion","edu")
  cur_data$sex=as.factor(cur_data$sex)
  cur_data$scanner=as.factor(cur_data$scanner)
  if (test_flip_sign[i] == 1) {
    cur_data$cog=-cur_data$cog
  }
  
  k_factor=3
  gm=gam(cog~s(s1,k=k_factor)+s(f1,k=k_factor)+s(s2,k=k_factor)+s(f2,k=k_factor)+s(s3,k=k_factor)+s(f3,k=k_factor)+age+sex+edu+scanner+motion,data=cur_data)
  summary(gm)
  
  all_models[[i]]=gm
  s=summary(gm)
  model_devs=rbind(model_devs,s$dev.expl)
  term_ts=rbind(term_ts,s$p.table[,3])
  term_t_ps=rbind(term_t_ps,s$p.table[,4])
  term_fs=rbind(term_fs,s$s.table[,3])
  term_f_ps=rbind(term_f_ps,s$s.table[,4])
  term_f_dfs=rbind(term_f_dfs,s$s.table[,1])
  all_preds[[i]]=gm$fitted.values
}

#write.table(term_fs,file="gam_fs.txt",sep=" ",col.names = FALSE,row.names=FALSE)
#write.table(term_f_ps,file="gam_f_ps.txt",sep=" ",col.names = FALSE,row.names=FALSE)
