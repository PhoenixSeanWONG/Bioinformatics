## Relationship between PFS/OS and other variable（gender, subtypes, tumour_purity, disease_status, anat_loc, Treatment, Response, AvgSpotLen, RECIST）

```{r}
# Convert Dead/Alive(Dead = True) to vital_status and to a numeric value (Alive = 1, Dead = 0)
colnames(dead_patients) <- c("Patient_ID", "vital_status", "PMID")

dead_patients <- dead_patients %>%
  mutate(vital_status = 0)
```

```{r}
colnames(alive_patients) <- c("Patient_ID", "vital_status", "PMID")

alive_patients <- alive_patients %>%
  mutate(vital_status = 1)
```

```{r}
updated_vital_status <- bind_rows(dead_patients, alive_patients)

final_combined_withHugo <- final_combined_withHugo %>%
  select(-vital_status)

final_combined_withHugo <- final_combined_withHugo %>%
  left_join(updated_vital_status %>% select(Patient_ID, vital_status), by = "Patient_ID")

table(final_combined_withHugo$vital_status)
```

### Check data format and modify
```{r}
str(final_combined_withHugo)
final_combined_withHugo$vital_status <- as.numeric(final_combined_withHugo$vital_status)
```

```{r}
names(final_combined_withHugo)[names(final_combined_withHugo) == "Response (0/1)"] <- "Response"
```

```{r}
str(final_combined_withHugo)
# Convert all specified variables to factor types
factor_vars <- c("gender", "PFS_Group", "OS_Group", "subtypes", "tumour_purity", "disease_status", "anat_loc", "Treatment", "Response", "AvgSpotLen", "RECIST")
final_combined_withHugo[factor_vars] <- lapply(final_combined_withHugo[factor_vars], as.factor)
```

#### Group
```{r}
# Convert Group variables to factor types
final_combined_withHugo$PFS_Group <- as.factor(final_combined_withHugo$PFS_Group)

final_combined_withHugo$OS_Group <- as.factor(final_combined_withHugo$OS_Group)

# Check number of levels in Group
table(final_combined_withHugo$PFS_Group)
table(final_combined_withHugo$OS_Group)
```

#### subtypes
```{r}
# Convert subtypes variables to factor types
final_combined_withHugo$subtypes <- as.factor(final_combined_withHugo$subtypes)
# Check number of levels in subtypes
table(final_combined_withHugo$subtypes)
```

#### tumour_purity
```{r}
# Convert tumour_purity variable to factor type
final_combined_withHugo$tumour_purity <- as.factor(final_combined_withHugo$tumour_purity)
# Check number of levels in tumour_purity
table(final_combined_withHugo$tumour_purity)
```

#### disease_status
```{r}
# Convert disease_status variable to factor type
final_combined_withHugo$disease_status <- as.factor(final_combined_withHugo$disease_status)
# Check number of levels in disease_status
table(final_combined_withHugo$disease_status)
```

#### anat_loc
```{r}
# Convert anat_loc variable to factor type
final_combined_withHugo$anat_loc <- as.factor(final_combined_withHugo$anat_loc)
# Check level number in anat_loc
table(final_combined_withHugo$anat_loc)
```

#### Treatment
```{r}
# Convert Treatment to lowercase
final_combined_withHugo$Treatment <- tolower(final_combined_withHugo$Treatment)

# Convert Treatment variables to factor types
final_combined_withHugo$Treatment <- as.factor(final_combined_withHugo$Treatment)

# Check number of levels in Treatment
table(final_combined_withHugo$Treatment)
```

#### Response
```{r}
# Convert Response variables to factor types
final_combined_withHugo$Response <- as.factor(final_combined_withHugo$Response)
# Check number of levels in Response
table(final_combined_withHugo$Response)
```

#### AvgSpotLen
```{r}
# Convert AvgSpotLen variables to factor types
final_combined_withHugo$AvgSpotLen <- as.factor(final_combined_withHugo$AvgSpotLen)
# Check number of levels in AvgSpotLen
table(final_combined_withHugo$AvgSpotLen)
```

#### RECIST
```{r}
# Convert RECIST variables to factor types
final_combined_withHugo$RECIST <- as.factor(final_combined_withHugo$RECIST)
# Check number of levels in RECIST
table(final_combined_withHugo$RECIST)
```

```{r}
str(final_combined_withHugo)
```

```{r}
# Check for missing values
colSums(is.na(final_combined_withHugo))
```
