# Bacillus-Pseudomonas biofilm evolution population sequencing analysis
# 1. Install and load the packages, then read the mutation table
packages <- c("tidyverse", "readxl","openxlsx","ggplot2","ggpubr","ggsci","vegan","ape","gridExtra","cowplot","multcompView","agricolae","reshape2","lmerTest","brms","emmeans")
install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)

# Read the mutation table
df <- read_excel("Mutation.all.xlsx",sheet = 1,col_names = T)


# -------------------------------------------------------------------------------------
# 2. Filter the mutation table based on the following criteria
# (1) Filter out the mutations found in the ancestors
# (2) Filter out the mutations that never reach a cumulative frequency of 5%
# (3) Filter out synonymous mutations
df.filter <- df%>%
  filter(T0==0)%>%
  filter(max>=0.05)%>%
  filter(snp_type!="synonymous")

# The filtered mutations include fixed, transient, and fluctuating mutations.
df1 <- df.filter%>%
  select(Lineage,Group,T0,T1,T3,T6,T9,T12,T15,T18,position,Class,snp_type)%>%
  pivot_longer(names_to = "Transfer",  # Convert column names (T0, T1, ...) into a "Transfer" column
               values_to = "Freq",     # Convert the values into the "Freq" column
               cols = starts_with("T"))   #Apply to all columns starting with "T"

# Summarize mutations by lineage, group, and transfer time
mutsum <- df1 %>%
  filter(Freq != 0) %>%   # Only include rows with non-zero frequency
  filter(Group != "P-co") %>%   # Exclude "P-co" group            
  count(Lineage, Group, Transfer, position) %>%  # Count mutations per position
  mutate(Lineage = gsub("[BP]", "", Lineage))    # Remove "B" and "P" from Lineage labels

# Calculate the total number of mutations for each lineage by transfer, and plot as a line chart
lineage <- mutsum %>%
  group_by(Lineage, Group, Transfer) %>%
  summarise(total = sum(n)) %>%
  ungroup() %>%
  mutate(Transfer = factor(gsub("[T]", "", Transfer),   # Convert Lineage to a factor
                           levels = c("1", "3", "6", "9", "12", "15", "18"))) 

ggplot(lineage,aes(x = Transfer, y = total,group=Lineage,color=Lineage))+
  geom_line()+
  geom_point()+
  facet_wrap(~Group)+
  labs(x = "Transfer",y = "Average number of cumulative\nmutations across transfers")+
  theme_bw()+
  theme(axis.text = element_text(size=6,color="black"),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.position = "right",
        legend.title = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key.size = unit(0.9, "lines"),legend.margin = margin(0, 0, 0, 0))
ggsave("Number_of_mutations (all filter).pdf",width = 85,height = 50,units = "mm")

# Calculate mutation types for each lineage and plot as a bar chart
ave1 <- df1%>%
  filter(Freq != 0) %>%
  count(Lineage, Group, Transfer, snp_type) %>%    # Count mutations by snp_type
  group_by(Lineage, Group, snp_type) %>%
  summarise(ave = mean(n), .groups = "drop") %>%  # Calculate the average number of mutations
  filter(Group != "P-co")                         # Exclude "P-co" group

ggplot(ave1, aes(x = Lineage, y = ave, fill = snp_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Lineage",
       y = "Average number of cumulative\nmutations across transfers",
       fill = "Mutation Type") +
  facet_wrap(~Group, scales = "free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(size=6,color="black",angle = 90),
        axis.text.y = element_text(size=6,color="black"),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.position = "right",
        legend.title = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key.size = unit(0.9, "lines"),legend.margin = margin(0, 0, 0, 0)) 
ggsave("Mutation_Type (all filter).pdf",width = 95,height = 50,units = "mm")

 
# Calculate Jaccard similarity index
# Function to calculate Jaccard index
calculate_jaccard <- function(df) {
  combinations <- combn(unique(df$Lineage), 2, simplify = FALSE)
  jaccard_indices <- lapply(combinations, function(pair) {
    lineage1 <- pair[1]
    lineage2 <- pair[2]
    positions_lineage1 <- df %>%
      filter(Lineage == lineage1) %>%
      distinct(position)
    positions_lineage2 <- df %>%
      filter(Lineage == lineage2) %>%
      distinct(position)
    intersect_length <- length(intersect(positions_lineage1$position, positions_lineage2$position))
    union_length <- length(union(positions_lineage1$position, positions_lineage2$position))
    jaccard_index <- intersect_length / union_length
    return(data.frame(Lineage1 = lineage1, Lineage2 = lineage2, JaccardIndex = jaccard_index))
  })
  return(bind_rows(jaccard_indices))
}

# Calculate Jaccard index for each group
jaccard_results1 <- df1 %>%
  filter(Group!="P-co",Transfer!="T0")%>%
  group_by(Group) %>%
  nest() %>%
  mutate(jaccard_data = map(data, calculate_jaccard)) %>%
  select(-data) %>%
  unnest(jaccard_data)

# Plot using ggplot2 with stat_compare_means
ggplot(jaccard_results1, aes(x = Group, y = JaccardIndex, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  labs(title = "Jaccard Indices for pairwise Lineages ",
       x = "Group",
       y = "Jaccard similarity of fixed, transient,\nand fluctuating mutations") +
  stat_compare_means(method = "t.test", comparisons = list(c("B-co", "B-mono")), label = "p.signif") +
  ylim(0,0.4)+
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black"),
        axis.title = element_text(size=7),
        legend.position = "none",
        plot.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("Jaccard_Index (all filter).pdf",width = 40,height =50,units="mm")
ggsave("Jaccard_Index (all filter).png",width = 40,height =50,units="mm")




# ---------------------------------------------------------------------------------
# 3.The above plots showed lineage variability, which may result from transient mutations which occurred in a single transfer time
# In the following analysis, filter out the transient mutations.
df2 <- df.filter%>%
  filter(occurrence!="1")%>%     # Filter out transient mutations
  select(Lineage,Group,T0,T1,T3,T6,T9,T12,T15,T18,position,Class,snp_type)%>%
  pivot_longer(names_to = "Transfer", values_to = "Freq", cols = starts_with("T"))

jaccard_results2 <- df2 %>%
  filter(Group!="P-co",Transfer!="T0")%>%
  group_by(Group) %>%
  nest() %>%
  mutate(jaccard_data = map(data, calculate_jaccard)) %>%
  select(-data) %>%
  unnest(jaccard_data)
ggplot(jaccard_results2, aes(x = Group, y = JaccardIndex, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  labs(title = "Jaccard Indices for pairwise Lineages ",
       x = "Group",
       y = "Jaccard similarity of\nnon-transient mutations") +
  stat_compare_means(method = "t.test", comparisons = list(c("B-co", "B-mono")), 
                     label = "p.signif",label.y=0.56,size=3) +
  ylim(0.15,0.6)+
  theme_bw() +
  theme(axis.text = element_text(size=7, color="black"),
        axis.title = element_text(size=7),
        legend.position = "none",
        plot.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("Jaccard_Index (non-transient).pdf",width = 40,height =50,units="mm")
ggsave("Jaccard_Index (non-transient).png",width = 40,height =50,units="mm")




# -------------------------------------------------------------------------------------
# 4. Plot the accumulation of non-transient mutations, and do statistic analysis
mutsum2 <- df2 %>%
  filter(Freq != 0) %>%
  filter(Group != "P-co")%>%
  count(Lineage, Group, Transfer, position) %>% 
  mutate(Lineage = gsub("[BP]", "", Lineage))

lineage2 <- mutsum2 %>%
  group_by(Lineage, Group, Transfer) %>%
  summarise(total = sum(n)) %>%
  ungroup()  %>%
  mutate(Transfer = factor(gsub("[T]", "", Transfer), 
                           levels = c("1", "3", "6", "9", "12", "15", "18"))) 

ggplot(lineage2,aes(x = Transfer, y = total,group=Lineage,color=Lineage))+
  geom_line()+
  geom_point()+
  facet_wrap(~Group)+
  labs(x = "Transfer",y = "Cumulative number of\nnon-transient mutations")+
  theme_bw()+
  theme(axis.text = element_text(size=6,color="black"),
        axis.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.position = "right",
        legend.title = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key.size = unit(0.9, "lines"),legend.margin = margin(0, 0, 0, 0))+
  ylim(0,35)
ggsave("Number_of_mutations (non transient).pdf",width = 90,height = 50,units = "mm")
ggsave("Number_of_mutations (non transient).png",width = 90,height = 50,units = "mm")


# Statistic analysis
# linear mixed model
model <- lmer(total ~ Group * Transfer + (1 | Lineage), data = lineage2)
summary(model)
anova_result <- anova(model)
confint(model, level = 0.95)
pairwise_comparisons <- emmeans(model, pairwise ~ Group | Transfer)
summary(pairwise_comparisons)
summary(model)$coefficients["GroupB-mono", ]
VarCorr(model)






# -------------------------------------------------------------------------------------
# 5. Calculate the average number of non-transient mutations of each lineage across transfers
# and plot the mutation type
ave2 <- df2%>%
  filter(Freq != 0) %>%
  count(Lineage, Group, Transfer, snp_type) %>%
  group_by(Lineage, Group, snp_type) %>%
  summarise(ave = mean(n), .groups = "drop") %>%
  filter(Group != "P-co")

ggplot(ave2, aes(x = Lineage, y = ave, fill = snp_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Lineage",
       y = "Average number of non-transient\nmutations across transfers",
       fill = "Mutation Type") +
  facet_wrap(~Group, scales = "free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(size=6,color="black",angle = 90),
        axis.text.y = element_text(size=6,color="black"),
        axis.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.position = "right",
        legend.title = element_text(size=7),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key.size = unit(0.9, "lines"),legend.margin = margin(0, 0, 0, -5)) +
  ylim(0,25)
ggsave("Mutation_Type(non transient).png",width = 70,height = 50,units = "mm")
ggsave("Mutation_Type(non transient).pdf",width = 70,height = 50,units = "mm")





# ---------------------------------------------------------------------------------
# 6. Calculate bray distance, and plot NMDS plot
df2_2 <- df2%>%
  group_by(Lineage,Group,Transfer,position)%>%
  summarise(Freqsum = sum(Freq))

# Define the conditions and their corresponding stages
conditions <- list(
  Early = c("T0", "T1", "T3"),
  Middle = c("T6","T9", "T12"),
  Late = c("T15", "T18"))

df3 <- df2_2 %>%
  filter(Group!="P-co",Transfer!="T0")%>%
  mutate(Stage = case_when(
    Transfer %in% conditions$Early ~ "Early",
    Transfer %in% conditions$Middle ~ "Middle",
    Transfer %in% conditions$Late ~ "Late",
    TRUE ~ NA_character_ ))%>%         # Handle other cases, if any
  pivot_wider(names_from = position,
              values_from = Freqsum,
              values_fill = 0)

# Select only the numeric columns
df3_numeric <- df3[, -(1:4)]
df3_numeric[] <- lapply(df3_numeric, function(x) as.numeric(as.character(x)))

# Calculate the Bray-Curtis distance matrix
dist_matrix <- vegdist(df3_numeric, method = "bray")

# Perform NMDS
nmds <- metaMDS(dist_matrix)

# Extract NMDS coordinates
nmds_coords <- as.data.frame(scores(nmds))

#Add Group column to NMDS coordinates
nmds_coords <- cbind(df3[, 1:4], nmds_coords)
nmds_coords <- nmds_coords%>%
  unite("GroupStage",c("Group","Stage"),sep = "_", remove = F)

# Extract stress value
stress_value <- nmds$stress
stress_value

# Extract eigenvalues
eigenvalues <- nmds$points

# Calculate the proportion of variance explained by each NMDS axis
var_explained <- eigenvalues^2 / sum(eigenvalues^2)

# Calculate the cumulative variance explained by NMDS1 and NMDS2
cumulative_var_explained <- cumsum(var_explained)

# Print the variance explained by NMDS1 and NMDS2
cat("Variance explained by NMDS1:", var_explained[1] * 100, "%\n")
cat("Variance explained by NMDS2:", var_explained[2] * 100, "%\n")
cat("Cumulative variance explained by NMDS1 and NMDS2:", cumulative_var_explained[2] * 100, "%\n")

# Perform ANOVA for NMDS1
anova_result1 <- aov(NMDS1 ~ GroupStage, data = nmds_coords)
summary(anova_result1)

# Perform Tukey's HSD posthoc test
tukey_result1 <- HSD.test(anova_result1, "GroupStage", group = TRUE)

# Extract the compact letter display
cld1 <- tukey_result1$groups
cld1 <- as.data.frame(cld1)
cld1$GroupStage <- rownames(cld1)

# Merge CLD with the original data frame for plotting
nmds_coords <- nmds_coords %>%
  left_join(cld1, by = c("GroupStage" = "GroupStage"))

# Perform ANOVA for NMDS2
anova_result2 <- aov(NMDS2 ~ GroupStage, data = nmds_coords)
summary(anova_result2)

# Perform Tukey's HSD posthoc test
tukey_result2 <- HSD.test(anova_result2, "GroupStage", group = TRUE)

# Extract the compact letter display
cld2 <- tukey_result2$groups
cld2 <- as.data.frame(cld2)
cld2$GroupStage <- rownames(cld2)

# Merge CLD with the original data frame for plotting
nmds_coords <- nmds_coords %>%
  left_join(cld2, by = c("GroupStage" = "GroupStage"))

# Plot NMDS using ggplot
nmds_coords$Transfer <- factor(nmds_coords$Transfer, levels = c("T0", "T1", "T3", "T6", "T9", "T12", "T15", "T18"))
nmds_coords$Group <- factor(nmds_coords$Group, levels = c( "B-co","B-mono"))
nmds_coords$Stage <- factor(nmds_coords$Stage, levels = c("Early", "Middle","Late"))
nmds_coords$GroupStage <- factor(nmds_coords$GroupStage, levels = c("B-mono_Early","B-mono_Middle","B-mono_Late","B-co_Early","B-co_Middle","B-co_Late")) 

ggplot(nmds_coords, aes(x = NMDS1.x, y = NMDS2.x, color = Group,shape=Stage)) +
  geom_point() +
  labs(x="NMDS1",y="NMDS2",title = "NMDS based on Bray-Curtis dissimilarity")+
  scale_shape_manual(values = c(1,9, 17)) +
  theme_bw() +
  theme(axis.text = element_text(size = 7, color = "black"),axis.title = element_text(size = 7),
        legend.text = element_text(size = 7),legend.position = "right",
        legend.title = element_text(size = 7),plot.title=element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "lines"),legend.margin = margin(0, 0, 0, -5))

ggsave("NMDS (non-transient).pdf",width = 63,height = 50,units = "mm")
ggsave("NMDS (non-transient).png",width = 63,height = 50,units = "mm")


# Perform PERMANOVA using the 'Group' column
permanova_result <- adonis2(dist_matrix ~ Group, data = df3, permutations = 999)
# Print the results
print(permanova_result)

# Extract and calculate the percentage of variance explained
r_squared <- permanova_result$R2[1]  # R2 for the grouping factor
percentage_variance_explained <- r_squared * 100
# Print the percentage of variance explained
cat("Percentage of variance explained by the grouping factor:", percentage_variance_explained, "%\n")




# ----------------------------------------------------------------------------------------
# 7. Plot the bray distance over time
# Convert the distance matrix to a data frame
dist_matrix_df <- as.data.frame(as.matrix(dist_matrix))

# Create a combined identifier for each sample
df3$SampleID <- paste(df3$Group, df3$Transfer, df3$Lineage, sep = "_")

# Add the row and column names to the distance matrix
rownames(dist_matrix_df) <- df3$SampleID
colnames(dist_matrix_df) <- df3$SampleID

# Melt the distance matrix to a long format
dist_long <- melt(as.matrix(dist_matrix_df), varnames = c("SampleID1", "SampleID2"), 
                  value.name = "Distance")%>% 
  filter(SampleID1 != SampleID2) %>%
  separate(SampleID1, into = c("Group1", "Transfer1", "Lineage1"), sep = "_") %>%
  separate(SampleID2, into = c("Group2", "Transfer2", "Lineage2"), sep = "_")


# Filter for comparisons
dist_compare<- dist_long %>% 
  filter(Group1 != Group2,Transfer1 == Transfer2,Group1 == "B-mono")%>%
  mutate(Transfer1 = factor(gsub("[T]", "", Transfer1), levels = c("1", "3", "6", "9", "12", "15", "18"))) %>%
  mutate(Transfer2 = factor(gsub("[T]", "", Transfer2), levels = c("1", "3", "6", "9", "12", "15", "18"))) 

modelcompare <- lm(Distance ~ Transfer1, data = dist_compare)
model_compare <- summary(modelcompare)
R_squared <- model_compare$r.squared
p_values <- summary(modelcompare)$coefficients[grepl("Transfer1", rownames(summary(modelcompare)$coefficients)), "Pr(>|t|)"]
R <- signif(sqrt(R_squared), digits = 2)


# Plot the distances
ggplot(dist_compare, aes(x = Transfer1, y = Distance)) +
  geom_jitter(aes(color = Transfer1), width = 0.15, size = 0.5) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  labs(x = "Transfer", y = "Bray-Curtis Distance\nBetween B-mono and B-co") +
  annotate("text", x = -Inf, y = Inf, 
           label = paste("Correlation coefficient: ", R, 
                         "\n  R²: ", signif(R_squared, digits = 2), 
                         "\n  p: ", signif(min(p_values), digits = 3)), 
           hjust = -0.05, vjust = 1.1, size = 2, color = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 6, color = "black"),
        axis.title = element_text(size = 7),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("BCdistance_compare (non-transient).pdf",units = "mm",width=48,height=50)
ggsave("BCdistance_compare (non-transient).png",units = "mm",width=48,height=50)



# Filter within B-mono group
dist_mono <- dist_long %>%
  filter(Transfer1 == Transfer2,Group1 == "B-mono",Group2 == "B-mono")%>%
  mutate(Transfer1 = factor(gsub("[T]", "", Transfer1), levels = c("1", "3", "6", "9", "12", "15", "18"))) %>%
  mutate(Transfer2 = factor(gsub("[T]", "", Transfer2), levels = c("1", "3", "6", "9", "12", "15", "18"))) 

modelmono <- lm(Distance ~ Transfer1, data = dist_mono)
model_mono <- summary(modelmono)
R_squared <- model_mono$r.squared
p_values <- summary(modelmono)$coefficients[grepl("Transfer1", rownames(summary(modelmono)$coefficients)), "Pr(>|t|)"]
R <- signif(sqrt(R_squared), digits = 2)

ggplot(dist_mono, aes(x = Transfer1, y = Distance)) +
  geom_jitter(aes(color = Transfer1), width = 0.15, size = 0.5) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  labs(x = "Transfer", y = "Bray-Curtis Distance\nwithin B-mono") +
  annotate("text", x = -Inf, y = Inf, 
           label = paste("Correlation coefficient: ", R, 
                         "\n  R²: ", signif(R_squared, digits = 2), 
                         "\n  p: ", signif(min(p_values), digits = 3)), 
           hjust = -0.05, vjust = 1.1, size = 2, color = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("BCdistance_mono (non-transient).pdf",units = "mm",width=60,height=55)


# Filter within B-co group
dist_co <- dist_long %>%
  filter(Transfer1 == Transfer2,Group1 == "B-co",Group2 == "B-co")%>%
  mutate(Transfer1 = factor(gsub("[T]", "", Transfer1), levels = c("1", "3", "6", "9", "12", "15", "18"))) %>%
  mutate(Transfer2 = factor(gsub("[T]", "", Transfer2), levels = c("1", "3", "6", "9", "12", "15", "18"))) 

modelco <- lm(Distance ~ Transfer1, data = dist_co)
model_co <- summary(modelco)
R_squared <- model_co$r.squared
p_values <- summary(modelco)$coefficients[grepl("Transfer1", rownames(summary(modelco)$coefficients)), "Pr(>|t|)"]
R <- signif(sqrt(R_squared), digits = 2)

ggplot(dist_co, aes(x = Transfer1, y = Distance)) +
  geom_jitter(aes(color = Transfer1), width = 0.15, size = 0.5) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "red", linewidth = 0.8) +
  labs(x = "Transfer", y = "Bray-Curtis Distance\nwithin B-co") +
  annotate("text", x = -Inf, y = Inf, 
           label = paste("Correlation coefficient: ", R, 
                         "\n  R²: ", signif(R_squared, digits = 2), 
                         "\n  p: ", signif(min(p_values), digits = 3)), 
           hjust = -0.05, vjust = 1.1, size = 2, color = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("BCdistance_co (non-transient).pdf",units = "mm",width=60,height=55)



# -------------------------------------------------------------------------------------
# 8. Filter the fixed mutations, and plot for each lineage
fixed <- df.filter%>%
  filter(classification!="unfixed")%>%     # Filter out unfixed mutations
  select(Group,Lineage,T0:T18,Gene,snp_type)%>%
  pivot_longer(names_to = "Transfer", 
               values_to = "Freq",     
               cols = starts_with("T"))%>%
  unite("Gene_Mutation",c("Gene","snp_type"),sep = "_", remove = F)%>%
  mutate(Transfer = factor(gsub("[T]", "", Transfer),   # Convert Lineage to a factor
                           levels = c("0","1", "3", "6", "9", "12", "15", "18"))) 

fixedmono <- fixed%>%filter(Group=="B-mono")
fixedco <- fixed%>%filter(Group=="B-co")
fixedpse <- fixed%>%filter(Group=="P-co")

# Function to create and save the plot
create_plot <- function(data, title, filename) {
  ggplot(data, aes(x = Transfer, y = Freq, group = Gene_Mutation, color = Gene_Mutation, linetype = snp_type)) +
    geom_line() +
    facet_wrap(~Lineage, ncol = 4) +
    labs(y = "Frequency", title = title) +
    scale_linetype_manual(values = c("intergenic" = "dotted", "InDel" = "solid", "missence" = "dashed")) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.position = "right",
      legend.title = element_text(size = 8),
      plot.title = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.9, "lines"),
      legend.margin = margin(0, 0, 0, 0)
    )
  ggsave(filename, height = 9, width = 18, units = "cm")
}

# Create and save the plots
create_plot(fixedmono, "Mono-evolved Bacillus", "Fixed_Mono.pdf")
create_plot(fixedco, "Co-evolved Bacillus", "Fixed_Co.pdf")
create_plot(fixedpse, "Co-evolved Pseudomonas", "Fixed_Pse.pdf")
create_plot(fixedmono, "Mono-evolved Bacillus", "Fixed_Mono.png")
create_plot(fixedco, "Co-evolved Bacillus", "Fixed_Co.png")
create_plot(fixedpse, "Co-evolved Pseudomonas", "Fixed_Pse.png")





# ----------------------------------------------------------------------------------------
# 9. Filter out mutations occurred in intergenic regions, tRNA, rRNA
# Export as 24 excel files for each population (8 mono-evolved Bacillus, 8 co-evolved Bacilus, 8 co-evolved Pseudomonas), these are used as the input files of lolipop
dfmuller <- df.filter%>%
  filter(snp_type!="intergenic")%>%
  filter(snp_type!="tRNA")%>%
  filter(snp_type!="rRNA")%>%
  filter(occurrence >=2)

# Use gsub to remove the 'T' from column names that start with 'T'
names(dfmuller) <- gsub("^T", "", names(dfmuller))
  
lineages <- paste0("B", 1:8)
groups <- c("", "P")  # "" for B-mono and "P" for BP-mono and Pseudomonas

# Loop through each group and lineage
for (group in groups) {
  for (lineage in lineages) {
    # Create lineage identifier
    lineage_id <- paste0(group, lineage)
    
    # Filter df based on lineage and chromosome condition
    if (group == "") {
      # Bacillus monoculture populations
      data_to_save <- subset(dfmuller, Lineage == lineage_id)
    } else if (group == "P") {
      # Bacillus co-culture with Pseudomonas populations
      data_to_save <- subset(dfmuller, Lineage == lineage_id, Chromosome == "CP006890")
      write.xlsx(data_to_save, paste0(lineage_id, ".xlsx"))
      
      # Pseudomonas populations
      data_to_save <- subset(dfmuller, Lineage == lineage_id, Chromosome == "CP046538")
      write.xlsx(data_to_save, paste0("P", lineage, ".xlsx"))
      next
    }
    
    # Save the filtered data to an Excel file
    write.xlsx(data_to_save, paste0(lineage_id, ".xlsx"))
  }
}

# These files are used as input file for lolipop
# lolipop requires linux system
# We used the default parameter to infer genotype and genealogy of each population, e.g.:
# lolipop lineage --input dataset/B1.xlsx  --output B1






# ---------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------
# Not shown in the figures
# Shannon diversity (richness+eveness)
# Calculate Shannon diversity
shannon_diversity <- df2 %>%
  group_by(Transfer, Group, Lineage) %>%
  summarize(Shannon_Diversity = diversity(Freq, index = "shannon"))
shannon_diversity$Transfer <- factor(shannon_diversity$Transfer,levels=c("T0","T1","T3","T6","T9","T12","T15","T18"))
# Plot using ggplot with facet_wrap
ggplot(shannon_diversity, aes(x = Transfer, y = Shannon_Diversity, color = Lineage, group = Lineage)) +
  geom_line() +
  geom_point() +
  labs(x = "Transfer", y = "Shannon diversity of mutations") +
  facet_wrap(~Group)+
  theme_bw()+
  theme(axis.text = element_text(size=8,color="black"),
        axis.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.position = "right",
        legend.title = element_text(size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("Shannon_index (non-transient).png",width = 150,height = 70,units = "mm")



