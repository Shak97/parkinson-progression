library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(catboost)
library(caret)       # for general model fitting
library(rpart)       # for fitting decision trees
library(ipred)       # for fitting bagged decision trees



smape <- function(real, pred) {
  return(100/length(real) * sum(2 * abs(pred - real) / (abs(real) + abs(pred))))
}

directory <- 'C:\\course-work\\semester-1\\digital-content-and-data-science\\project\\parkinsons-disease-progression\\'
protein_file <- 'train_proteins.csv'
peptides_file <- 'train_peptides.csv'
clinical_file <- 'train_clinical_data.csv'

train_proteins <- read.csv(file.path(directory, protein_file))
train_peptides <- read.csv(file.path(directory, peptides_file))
train_clinical <- read.csv(file.path(directory, clinical_file))

# Plotting Clinical Data --------------------------------------------------

str(train_clinical)
# Filter data for patient with ID 1517
plot_df <- train_clinical[train_clinical$patient_id == 1517, ]
# Create the plot
ggplot(plot_df, aes(x = visit_month)) +
  geom_line(aes(y = updrs_1, color = "updrs_1")) +
  geom_line(aes(y = updrs_2, color = "updrs_2")) +
  geom_line(aes(y = updrs_3, color = "updrs_3")) +
  geom_line(aes(y = updrs_4, color = "updrs_4")) +
  labs(x = "visit_month", y = "UPDRS", title = "Patient 1517") +
  scale_color_manual(values = c("blue", "red", "green", "orange")) +
  theme_minimal()

# Plotting Protein Data ---------------------------------------------------
# We will pick a random patient_id and plot the first 40 Protein entries(UniProt)
# of the patient and their NPX value against patient's visit month(visit_month).

# Filter data for patient with ID 1517
pro_plot_df <- train_proteins[train_proteins$patient_id == 1517, ]

# Get unique protein list
protein_list <- unique(pro_plot_df$UniProt)
protein_list <- protein_list[1:20]

# Filter pro_plot_df based on protein_list
pro_plot_df <- pro_plot_df[pro_plot_df$UniProt %in% protein_list, ]

# Get unique months and sort
unique_month <- unique(pro_plot_df$visit_month)
unique_month <- sort(unique_month)

# Calculate the number of rows for subplots
n_rows <- ceiling(length(unique_month) / 2)

# Create a list of plots
plot_list <- list()

# Iterate over unique_month
for (i in seq_along(unique_month)) {
  month <- unique_month[i]
  plot_df <- pro_plot_df[pro_plot_df$visit_month == month, ]
  p <- ggplot(plot_df, aes(x = UniProt, y = NPX)) +
    geom_bar(stat = "identity", fill = "blue") +
    coord_flip() +
    xlab("UniProt") +
    ylab("NPX") +
    ggtitle(paste("visit_month", month)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
  
  plot_list[[i]] <- p
}

# Combine plots using patchwork
combined_plot <- patchwork::wrap_plots(plotlist = plot_list, ncol = 2)

print(combined_plot)


preProcessData <- function(train_proteins, train_peptides){
  # Group by 'visit_id' and 'UniProt', and calculate the mean of 'NPX'
  df_protein_grouped <- train_proteins %>%
    group_by(visit_id, UniProt) %>%
    summarise(NPX = mean(NPX)) %>%
    ungroup()
  
  # Group by 'visit_id' and 'Peptide', and calculate the mean of 'PeptideAbundance'
  df_peptide_grouped <- train_peptides %>%
    group_by(visit_id, Peptide) %>%
    summarise(PeptideAbundance = mean(PeptideAbundance)) %>%
    ungroup()
  
  df_protein <- dcast(df_protein_grouped, visit_id ~ UniProt, value.var = "NPX")
  colnames(df_protein) <- c(names(df_protein)[1], colnames(df_protein)[2:length(colnames(df_protein))])
  df_protein <- as.data.frame(df_protein)
  
  df_peptide <- dcast(df_peptide_grouped, visit_id ~ Peptide, value.var = "PeptideAbundance")
  colnames(df_peptide) <- c(names(df_peptide)[1], colnames(df_peptide)[2:length(colnames(df_peptide))])
  df_peptide <- as.data.frame(df_peptide)
  
  
  pp_df <- merge(df_protein, df_peptide, by = "visit_id", all.x = TRUE)
  
  return(pp_df) 
}

handleMissigValues <- function(df){
  # Fill missing values with column means
  for (col in colnames(df)) {
    df[is.na(df[, col]), col] <- mean(df[, col], na.rm = TRUE)
  }
  return(df)
}

# call preProcessData
# args:
#       train_protiens
#       train_peptides
pp_df <- preProcessData(train_proteins, train_peptides)

# Fill missing values with column means
for (col in colnames(pp_df)) {
  pp_df[is.na(pp_df[, col]), col] <- mean(pp_df[, col], na.rm = TRUE)
}

# Create necessary variables
FEATURES  <- names(pp_df)[!(names(pp_df) %in% c("visit_id"))]
FEATURES <- c(FEATURES, "visit_month")
target <- c("updrs_1", "updrs_2", "updrs_3", "updrs_4")

smape_s <- c()
smape_train_lst <- c()
for (label in target){
  
  dataset_df <- merge(pp_df, train_clinical[c('visit_id', 'patient_id', 'visit_month', label)], by = 'visit_id', all.x = TRUE)
  dataset_df <- dataset_df[!is.na(dataset_df[,c(label)]), ]
  
  feature_list <- FEATURES
  feature_list <- c(feature_list, label)
  current_dataset <- dataset_df[,feature_list]
  
  #writing cleaned data files
  write.csv(current_dataset, file.path(directory, paste('cleaned_data_', label, '.csv')))
  
  set.seed(123)
  train_ratio <- 0.8
  train_indices <- sample(nrow(current_dataset), floor(train_ratio * nrow(current_dataset)))
  train_df <- current_dataset[train_indices, ]
  test_df <- current_dataset[-train_indices, ]
  
  print('---------------preprocessing done---------------')
  
  #create the formula using as.formula and paste
  formula <- as.formula(paste(label, ' ~ .' ))
  
  bagged_regressor <- bagging(
    formula = formula,
    data = train_df,
    nbagg = 100,  
    coob = TRUE,
    control = rpart.control(minsplit = 2, cp = 0)
  )
  
  print('---------------training done---------------')
  
  predictions <- predict(bagged_regressor, newdata = test_df)
  predictions_train <- predict(bagged_regressor, newdata = train_df)
  
  print('---------------predictions done---------------')
  
  bagging_smape <- smape(test_df[, label], predictions)
  smape_train <- smape(train_df[,label], predictions_train)
  
  print('---------------bagging_smape done---------------')
  
  smape_s <- c(smape_s, bagging_smape)
  smape_train_lst <- c(smape_train_lst, smape_train)
  
  print('---------------bagging_smape appending done---------------')
  
  print(bagging_smape)
  print(smape_train)
}

print('--------------------------------------------------------')

print(smape_s)
print(smape_train_lst)
