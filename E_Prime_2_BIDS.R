#This script automates the extraction of onset times from E-Prime to the BIDS format

## Preparation
set.seed(1234)
pacman::p_load(tidyverse, rprime,naniar)

## Import Data
### Get the names of all files
df <- as_data_frame(dir("BIDS_E-Prime", pattern = "txt", recursive = T))
df <- separate(data = df, value, into = c("Subject", "Task", "File"), sep = "/")
df$Subject <- as_factor(df$Subject)
df$Task <- as_factor(df$Task)
df$File <- as_factor(df$File)

df_fmri <- df %>% filter(Task %in% c("03_PS", "04_AL1", "05_AL2"))
df_fmri_no_excel <- df %>% filter(Task %in% c("03_PS", "04_AL1", "05_AL2")) %>% filter(!str_detect(File, "excel.txt"))
filename <- str_c("BIDS_E-Prime",df_fmri_no_excel$Subject[1],df_fmri_no_excel$Task[1],df_fmri_no_excel$File[1], sep = "/")

experiment_lines <- read_eprime(filename)
experiment_data <- FrameList(experiment_lines)


df_fmri_excel <- df %>% filter(Task %in% c("03_PS", "04_AL1", "05_AL2")) %>% filter(str_detect(File, "excel.txt"))
PS_excel <- df_fmri_excel %>% filter(Task == "03_PS")
AL1_excel <- df_fmri_excel %>% filter(Task == "04_AL1")
AL2_excel <- df_fmri_excel %>% filter(Task == "05_AL2")

### path to directory
path <- "/Users/marcocaviezel/Dropbox/Machine_Learning/R_stuff/E_Prime/BIDS_E-Prime"

### Process pattern separation

for(s in 1:length(PS_excel[[2]])){

  tmp_PS <- read.table(paste(path, PS_excel[[s,1]], PS_excel[[s,2]], PS_excel[[s,3]], sep ="/"),fileEncoding = 'UTF-16', skip=1, sep='\t', stringsAsFactors = F, header = T)
  current_subject = paste(path, PS_excel[[s,1]], PS_excel[[s,2]], PS_excel[[s,3]], sep ="/")
  print(current_subject)
  

  onset_rt_tmp <- tmp_PS %>% select(CueName, Stimulus, Art, Procedure.Block.,Comparison.ACC, Comparison1.ACC, WarmUpScanner.RTTime, Comparison.OnsetTime, Comparison1.OnsetTime,
                                   Comparison.RESP, Comparison1.RESP,Comparison.RT, Comparison1.RT, Comparison.OnsetToOnsetTime, Comparison1.OnsetToOnsetTime, 
                                   DiffFix2ExpBegin, DiffCompareExpBegin) %>% 
  mutate_at(c("Comparison.OnsetTime", "Comparison1.OnsetTime"), funs(ifelse(is.na(.), "x", .))) %>% 
  mutate(beh_onset = ifelse(Comparison.OnsetTime == "x", Comparison1.OnsetTime, Comparison.OnsetTime)) %>% 
  mutate(onset = as.numeric(beh_onset)-as.numeric(WarmUpScanner.RTTime)) %>%
  mutate(trial_type = ifelse(Procedure.Block. == "procvis"& Art == 1,"Visual_Same",
                                                 ifelse(Procedure.Block. == "procvis"& Art == 2,"Visual_Lure",
                                                        ifelse(Procedure.Block. == "procvis"& Art == 3,"Visual_Different",
                                                               ifelse(Procedure.Block. == "procaud" & Art == 1, "Auditory_Same",
                                                                      ifelse(Procedure.Block. == "procaud" & Art == 2, "Auditory_Lure", "Auditory_Different")))))) %>%
  mutate_at(c("Comparison.ACC", "Comparison1.ACC"), funs(ifelse(is.na(.), 9, .))) %>% mutate(accuracy = ifelse(Comparison.ACC+Comparison1.ACC==10, 1, 0)) %>% 
  mutate_at(c("Comparison.RESP", "Comparison1.RESP"), funs(ifelse(is.na(.), 9, .))) %>% mutate(button = ifelse(Comparison.RESP+Comparison1.RESP==10, "left", ifelse(Comparison.RESP+Comparison1.RESP==12, "right", NA))) %>% 
  mutate_at(c("Comparison.RT", "Comparison1.RT"), funs(ifelse(is.na(.), "x", .))) %>% mutate(response_time = ifelse(Comparison.RT == "x", Comparison1.RT, Comparison.RT)) %>% 
  mutate(response_time = ifelse(response_time==0,NA,response_time)) %>% 
  mutate_at(c("Comparison.OnsetToOnsetTime", "Comparison1.OnsetToOnsetTime"), funs(ifelse(is.na(.), "x", .))) %>% mutate(duration = ifelse(Comparison.OnsetToOnsetTime == "x", Comparison1.OnsetToOnsetTime, Comparison.OnsetToOnsetTime)) %>% 
  mutate(time_to_offset_alternative_onset = as.numeric(onset) + as.numeric(duration)) %>% 
  mutate(time_to_RT_alternative_onset_or_duration = as.numeric(onset) + as.numeric(response_time)) %>%
  mutate_at(c("onset", "accuracy", "response_time", "duration", "time_to_offset_alternative_onset", "time_to_RT_alternative_onset_or_duration"), funs(as.numeric(.))) %>% 
  mutate_at(c("onset", "response_time", "duration", "time_to_offset_alternative_onset", "time_to_RT_alternative_onset_or_duration"), funs(./1000)) %>% 
  select(onset, duration, trial_type, accuracy, button, response_time, time_to_offset_alternative_onset, time_to_RT_alternative_onset_or_duration)
  
  write.table(onset_rt_tmp, file = paste(path, PS_excel[[s,1]], PS_excel[[s,2]], paste(paste("events", PS_excel[[s,1]], PS_excel[[s,2]], sep = "_"),"tsv", sep = "."), sep ="/"), sep = "\t", row.names = F, quote = F)
  
  manual <- onset_rt_tmp %>% arrange(trial_type, onset) %>% filter(!is.na(button)) %>% filter(accuracy==1) %>% select(onset, trial_type)
  write.table(manual, file = paste(path, PS_excel[[s,1]], PS_excel[[s,2]], paste(paste("manual_events", PS_excel[[s,1]], PS_excel[[s,2]], sep = "_"),"tsv", sep = "."), sep ="/"), sep = "\t", row.names = F, quote = F)
}

### Process association learning (AL1)
for(s in 1:length(AL1_excel[[2]])){
  
tmp_AL1 <- read.table(paste(path, AL1_excel[[s,1]], AL1_excel[[s,2]], AL1_excel[[s,3]], sep ="/"),fileEncoding = 'UTF-16', skip=1, sep='\t', stringsAsFactors = F, header = T)
current_subject = paste(path, AL1_excel[[s,1]], AL1_excel[[s,2]], AL1_excel[[s,3]], sep ="/")
print(current_subject)

tmp2_AL1 <- tmp_AL1 %>% mutate(button_nr = rowSums(.[,c("FNP1.RESP", "FNP2.RESP","FNP3.RESP", "FNP4.RESP", "FNP5.RESP", "FNP6.RESP","FNP7.RESP", "FNP8.RESP","FNP9.RESP", "FNP10.RESP", "FNP11.RESP", "FNP12.RESP")],na.rm = T)) %>% 
  mutate(button = ifelse(button_nr == 3, "left", ifelse(button_nr == 1, "right",NA))) %>%
  mutate(scanner_start = (WarmUpScanner2.RTTime-2500)/1000) %>% 
  mutate(onset = ((as.numeric(rowSums(.[,c("FNP1.OnsetTime", "FNP2.OnsetTime","FNP3.OnsetTime", "FNP4.OnsetTime", "FNP5.OnsetTime", "FNP6.OnsetTime","FNP7.OnsetTime", "FNP8.OnsetTime","FNP9.OnsetTime", "FNP10.OnsetTime",
                                           "FNP11.OnsetTime", "FNP12.OnsetTime")],na.rm = T)))/1000)-scanner_start) %>% 
  mutate(response_time = (as.numeric(rowSums(.[,c("FNP1.RT", "FNP2.RT","FNP3.RT", "FNP4.RT", "FNP5.RT", "FNP6.RT","FNP7.RT", "FNP8.RT","FNP9.RT", "FNP10.RT", "FNP11.RT", "FNP12.RT")],na.rm = T))/1000)) %>% 
  replace_with_na(replace = list(response_time= 0)) %>% 
  mutate(time_to_RT_alternative_onset_or_duration = onset+response_time) %>% 
  mutate(time_to_offset_alternative_onset = (as.numeric(rowSums(.[,c("FNP1.OffsetTime", "FNP2.OffsetTime","FNP3.OffsetTime", "FNP4.OffsetTime", "FNP5.OffsetTime", "FNP6.OffsetTime","FNP7.OffsetTime", "FNP8.OffsetTime","FNP9.OffsetTime", "FNP10.OffsetTime", "FNP11.OffsetTime", "FNP12.OffsetTime")],na.rm = T))/1000)-scanner_start) %>% 
  mutate(duration = (as.numeric(rowSums(.[,c("FNP1.OnsetToOnsetTime", "FNP2.OnsetToOnsetTime","FNP3.OnsetToOnsetTime", "FNP4.OnsetToOnsetTime", "FNP5.OnsetToOnsetTime", "FNP6.OnsetToOnsetTime","FNP7.OnsetToOnsetTime", "FNP8.OnsetToOnsetTime","FNP9.OnsetToOnsetTime", "FNP10.OnsetToOnsetTime", "FNP11.OnsetToOnsetTime", "FNP12.OnsetToOnsetTime")],na.rm = T))/1000)) %>% 
  mutate(trial_type = ifelse(startsWith(Face, "C")==T, "Contingent", "NonContingent")) %>% 
  select(onset, duration, trial_type, button, response_time, time_to_offset_alternative_onset, time_to_RT_alternative_onset_or_duration)

write.table(tmp2_AL1, file = paste(path, AL1_excel[[s,1]], AL1_excel[[s,2]], paste(paste("events", AL1_excel[[s,1]], AL1_excel[[s,2]], sep = "_"),"tsv", sep = "."), sep ="/"), sep = "\t", row.names = F, quote = F)

manual <- tmp2_AL1 %>% arrange(trial_type, onset) %>% filter(!is.na(button)) %>% select(onset, trial_type)

write.table(manual, file = paste(path, AL1_excel[[s,1]], AL1_excel[[s,2]], paste(paste("manual_events", AL1_excel[[s,1]], AL1_excel[[s,2]], sep = "_"),"tsv", sep = "."), sep ="/"), sep = "\t", row.names = F, quote = F)

}

### Process association recall (AL2)
for(s in 1:length(AL2_excel[[2]])){
  
tmp_AL2 <- read.table(paste(path, AL2_excel[[s,1]], AL2_excel[[s,2]], AL2_excel[[s,3]], sep ="/"),fileEncoding = 'UTF-16', skip=1, sep='\t', stringsAsFactors = F, header = T)
current_subject = paste(path, AL2_excel[[s,1]], AL2_excel[[s,2]], AL2_excel[[s,3]], sep ="/")
print(current_subject)

tmp2_AL2 <- tmp_AL2 %>% 
  mutate(button_nr = rowSums(.[,c("Recall1.RESP", "Recall2.RESP","Recall3.RESP", "Recall4.RESP", "Recall5.RESP", "Recall6.RESP","Recall7.RESP", "Recall8.RESP","Recall9.RESP", "Recall10.RESP", "Recall11.RESP", "Recall12.RESP",
                                                        "Recall13.RESP", "Recall14.RESP", "Recall15.RESP", "Recall16.RESP", "Recall17.RESP")],na.rm = T)) %>% 
  mutate(button = ifelse(button_nr == 3, "left", ifelse(button_nr == 1, "right",NA))) %>% 
  mutate(scanner_start = (WaitForScannerRecall.RTTime-2500)/1000) %>% 
  mutate(onset = ((as.numeric(rowSums(.[,c("Recall1.OnsetTime", "Recall2.OnsetTime","Recall3.OnsetTime", "Recall4.OnsetTime", "Recall5.OnsetTime", "Recall6.OnsetTime","Recall7.OnsetTime", "Recall8.OnsetTime","Recall9.OnsetTime", "Recall10.OnsetTime", "Recall11.OnsetTime", "Recall12.OnsetTime",
                                           "Recall13.OnsetTime", "Recall14.OnsetTime", "Recall15.OnsetTime", "Recall16.OnsetTime", "Recall17.OnsetTime")],na.rm = T)))/1000)-scanner_start) %>% 
  mutate(response_time = (as.numeric(rowSums(.[,c("Recall1.RT", "Recall2.RT","Recall3.RT", "Recall4.RT", "Recall5.RT", "Recall6.RT","Recall7.RT", "Recall8.RT","Recall9.RT", "Recall10.RT", "Recall11.RT", "Recall12.RT",
                                                  "Recall13.RT", "Recall14.RT", "Recall15.RT", "Recall16.RT", "Recall17.RT")],na.rm = T))/1000)) %>% 
  replace_with_na(replace = list(response_time= 0)) %>% 
  mutate(time_to_RT_alternative_onset_or_duration = onset+response_time) %>%
  mutate(time_to_offset_alternative_onset = (as.numeric(rowSums(.[,c("Recall1.OffsetTime", "Recall2.OffsetTime","Recall3.OffsetTime", "Recall4.OffsetTime", "Recall5.OffsetTime", "Recall6.OffsetTime","Recall7.OffsetTime", "Recall8.OffsetTime","Recall9.OffsetTime", "Recall10.OffsetTime", "Recall11.OffsetTime", "Recall12.OffsetTime",
                                                                     "Recall13.OffsetTime", "Recall14.OffsetTime", "Recall15.OffsetTime", "Recall16.OffsetTime", "Recall17.OffsetTime")],na.rm = T))/1000)-scanner_start) %>%
  mutate(duration = (as.numeric(rowSums(.[,c("Recall1.OnsetToOnsetTime", "Recall2.OnsetToOnsetTime","Recall3.OnsetToOnsetTime", "Recall4.OnsetToOnsetTime", "Recall5.OnsetToOnsetTime", "Recall6.OnsetToOnsetTime","Recall7.OnsetToOnsetTime", "Recall8.OnsetToOnsetTime","Recall9.OnsetToOnsetTime", "Recall10.OnsetToOnsetTime", "Recall11.OnsetToOnsetTime", "Recall12.OnsetToOnsetTime",
                                             "Recall13.OnsetToOnsetTime", "Recall14.OnsetToOnsetTime", "Recall15.OnsetToOnsetTime", "Recall16.OnsetToOnsetTime", "Recall17.OnsetToOnsetTime")],na.rm = T))/1000)) %>% 
  mutate(trial_type = recode(Art, BLRL = "recently_learned_faces_baseline", BLFF = "famous_faces_baseline", RLFC = "recently_learned_faces_correct_combination", RLFI = "recently_learned_faces_false_combination",
                             FFC = "famous_faces_correct_combination", FFI = "famous_faces_false_combination")) %>% 
  mutate(trial_type_less_granular = recode(Art, BLRL = "recently_learned_faces_baseline", BLFF = "famous_faces_baseline", RLFC = "recently_learned_faces", RLFI = "recently_learned_faces", FFC = "famous_faces", FFI = "famous_faces" )) %>% 
  mutate(accuracy = rowSums(.[,c("Recall1.ACC", "Recall2.ACC","Recall3.ACC", "Recall4.ACC", "Recall5.ACC", "Recall6.ACC","Recall7.ACC", "Recall8.ACC","Recall9.ACC", "Recall10.ACC", "Recall11.ACC", "Recall12.ACC",
                                 "Recall13.ACC", "Recall14.ACC", "Recall15.ACC", "Recall16.ACC", "Recall17.ACC")],na.rm = T)) %>% 
  select(onset, duration, trial_type, accuracy, button, response_time, trial_type_less_granular, time_to_offset_alternative_onset, time_to_RT_alternative_onset_or_duration)

write.table(tmp2_AL2, file = paste(path, AL2_excel[[s,1]], AL2_excel[[s,2]], paste(paste("events", AL2_excel[[s,1]], AL2_excel[[s,2]], sep = "_"),"tsv", sep = "."), sep ="/"), sep = "\t", row.names = F, quote = F)

manual <- tmp2_AL2 %>% arrange(trial_type, onset) %>% filter(!is.na(button)) %>% filter(accuracy==1) %>% select(onset, trial_type)

write.table(manual, file = paste(path, AL2_excel[[s,1]], AL2_excel[[s,2]], paste(paste("manual_events", AL2_excel[[s,1]], AL2_excel[[s,2]], sep = "_"),"tsv", sep = "."), sep ="/"), sep = "\t", row.names = F, quote = F)

}


