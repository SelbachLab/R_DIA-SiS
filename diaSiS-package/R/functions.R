#' Convert Precursor.Id from DIA-NN 2.0 output
#' 
#' Converting Precursor.Id and Channel columns from the new output file back to the old format so the functions run. Do not use as stand alone. convert_20 wraps around it
#' 
#' @param precursor_id Precursor.Id
#' @param channel Channel
#' @returns Converted Precursor.Id
#' @import stringr
#' @export 
convert_precursor_id <- function(precursor_id = "Precursor.Id", channel = "Channel") {
  # Define the mapping of channels to old SILAC labels
  channel_map <- c("H" = "SILAC-H", "L" = "SILAC-L", "M" = "SILAC-M")
  
  # Extract sequence before (SILAC) and the modified amino acid
  match <- str_match(precursor_id, "(.+)([A-Z])\\(SILAC\\)(\\d*)")
  
  if(!is.na(match[1])){
    sequence <- match[2]  # Part before modified residue
    modified_aa <- match[3]  # Modified residue
    number <- match[4]  # Number after modification (if present)
    
    # Construct old format
    old_format <- paste0(sequence, "(", channel_map[channel], "-", modified_aa, ")", number)
    return(old_format)
  }else{
    return(precursor_id) # Fallback if no match found
  }
}

#' Convert DIA-NN 2.0 output to format in lower versions
#' 
#' Converting Precursor.Id and Channel columns from the new output file back to the old format so the functions run. Do not use as stand alone. convert.2.0 wraps around it
#' 
#' @param file file path
#' @returns tibble with Precursor.Id formatted as in earlier DIA-NN versions and removed Channel column.
#' @examples 
#' data <- convert_20(path/to/file/report.parquet)
#' @import stringr
#' @import dplyr
#' @import arrow
#' @export 
convert_20 <- function(file){
  read_parquet(file) %>% 
    mutate(Precursor.Id = map2_chr(Precursor.Id, Channel, convert_precursor_id)) %>% 
    select(-Channel) -> out
  
  return(out)
}





#' Clean DIA-NN output
#' 
#' Removing wrong charges and contaminants from any DIANN output file
#' 
#' @param file Path to report.tsv or report.parquet file or the already read in file
#' @param contaminats String defining contaminant protein groups
#' @returns DIA-NN output cleaned for contaminants and precursors with charge 1
#' @examples
#' data <- clean_DIANN(filepath)
#' @import readr
#' @import dplyr
#' @import arrow
#' @import stringr
#' @export 
clean_DIANN <- function(file, contaminants = "CON"){  
  if(str_detect(class(file), "character")){
    if(str_detect(file, ".tsv")){
      read_delim(file, delim = "\t") %>%
        filter(Precursor.Charge > 1, !str_detect(Protein.Group, contaminants)) %>%
        # create a precursor ID that only contains sequence and charge but not SILAC state
        mutate(Stripped.Sequence.Charge = str_remove_all(Precursor.Id, "\\(SILAC-[KR]-[LHM]\\)")) -> out
    }else{
      read_parquet(file) %>%
        filter(Precursor.Charge > 1, !str_detect(Protein.Group, contaminants)) %>%
        # create a precursor ID that only contains sequence and charge but not SILAC state
        mutate(Stripped.Sequence.Charge = str_remove_all(Precursor.Id, "\\(SILAC-[KR]-[LHM]\\)")) -> out
    }
  }else{
    file %>% 
      as_tibble() %>% 
      filter(Precursor.Charge > 1, !str_detect(Protein.Group, contaminants)) %>%
      # create a precursor ID that only contains sequence and charge but not SILAC state
      mutate(Stripped.Sequence.Charge = str_remove_all(Precursor.Id, "\\(SILAC-[KR]-[LHM]\\)")) -> out
  }

  return(out)
}



#' Flagging proteins, precursors and channels that pass filters
#' 
#' Creates a tibble with precursors and channels that pass certain filters. Filters used are the standard DIA-NN filters used in other publications (Global.PG.Q.Value < 0.01) and Channel.Q.Value < 0.03. 
#' 
#' @param data Output from clean_DIANN function or raw DIANN output. 
#' @param referenceChannel String defining which channel was spiked-in to use as a reference. Options are "H", "M" and "L". Default is H.
#' @param numberChannels Number of SILAC channels, default is 2
#' @param MBR Whether MBR was turned on and off. Default is off. Should only be turned on for Label free data
#' @param CalcCols Whether or not Ms1.Translated and Precursor.Translated must both have valid values or not. Default is "both" (non valid values are not left). This is recommended for eg timsTOF pro2 data analysed with DIA-NN 1.8.1. For sc data (and higher DIA-NN Version 1.8.2 beta 22 or higher) set to "any" (at least one of them have to be valid / uses Precursor.Normalised for higher DIA-NN versions). "Precursor.Translated" and "Ms1.Translated" can also be used seperately but is not recommenended
#' @param ChannelFilt Q value filter for Channel.Q.Value. If left NA, the default and suggested for DIA-NN 1.8.1 is 0.03. For higher versions the recommendation is ChannelFilt = 0.01
#' @param version DIA-NN version used to analyse the raw data. Default and recommended is 1.8.1. Set to 1.9 for higher versions. Recommended DIA-NN settings for 1.8.2 beta 22 or higher are QuantUMS: legacy; turn off MBR; and --channel-spec-norm
#' 
#' @returns Set of Protein.Groups and precursors per Run that passed filters based on an approach similar to requantify in Maxquant. If the reference passes the filter, other channels are also kept. Proteins identified with Requantify may not be quantified accurately. Therefore, the output also contains informations about how many channels per precursor passed the standard filter and an advanced strict filter that further reduces untrustworthy quantifications. Output contains: Run, Protein.Group, Stripped.Sequence.Charge (Precursor.Id without SILAC Channel information), Reference.Channel (Channel used for Requantify), Channels.Passed.Standard.Filter (Listing all Channels that passed the standard filter, not only the reference one), Channels.Passed.Strict.Filter (Listing all Channels that passed the strict filter, not only the reference one), N.Total.Precursors.Per.Protein.Standard (Number of requantify precursors that passed the standard filter), N.Precursors.With.Complete.SILAC.Channels.Per.Protein.Standard (Number of precursors per proteins where all channels passed the standard filter) and N.Precursors.With.Complete.SILAC.Channels.Per.Protein.Strict (Number of Precursors per Protein where all channels passed the strict filter).Function can also be used for non-SILAC runs. Set reference-Channel = "L" and numberChannels = 1.
#' @examples
#' data <- clean_DIANN(filepath)
#' 
#' # for bulk SILAC spike-in data
#' filter_frame <- filter_DIANN(data)
#' 
#'  # for sc SILAC spike-in data
#' filter_frame <- filter_DIANN(data, CalcCols = "any")
#' 
#'  # for DIA-NN 1.8.2 beta 22 or higher
#'  filter_frame <- filter_DIANN(data, CalcCols = "any", ChannelFilt = 0.01, version = "1.9")
#' @import readr
#' @import dplyr
#' @export filter_DIANN
#' 
filter_DIANN <- function(data, referenceChannel = "H", numberChannels = 2, MBR = "off", CalcCols = "both", ChannelFilt = NA, version = "1.8.1") {
  
  # Define threshold and default for ChannelFilt
  ChannelFilt <- ifelse(is.na(ChannelFilt), 0.03, ChannelFilt)
  filterChannelTh <- ifelse(numberChannels == 1, 0.01, ChannelFilt)
  
  # Set filter columns based on number of channels and MBR status
  filterChannel <- ifelse(numberChannels == 1, ifelse(MBR == "off", "Global.Q.Value", "Lib.Q.Value"), "Channel.Q.Value")
  filterPGglobal <- ifelse(MBR == "off", "Global.PG.Q.Value", "Lib.PG.Q.Value")
  searchChannel <- ifelse(numberChannels == 1, ".", "-[LMH]")
  
  # Set quantification columns based on CalcCols and version
  useCols <- if (CalcCols == "both" || CalcCols == "any") {
    if (version == "1.8.1") c("Ms1.Translated", "Precursor.Translated") else "Precursor.Normalised"
  } else {
    CalcCols
  }
  
  # Create precursor filter list with standardized filtering
  if (referenceChannel %in% c("H", "L", "M")) {
    
    # Set reference channel pattern based on number of channels
    referenceChannelPattern <- if (numberChannels == 1) ".*" else paste0("-", referenceChannel)
    
    # Filter data based on conditions
    data %>%
      filter(
        str_detect(Modified.Sequence, referenceChannelPattern),
        (!!sym(filterPGglobal)) < 0.01,
        (!!sym(filterChannel)) < filterChannelTh
      ) %>%
      filter_at(vars(useCols), if (CalcCols != "any") all_vars(!is.na(.) & . > 0) else any_vars(!is.na(.) & . > 0)) %>%
      group_by(Run, Protein.Group) %>%
      mutate(Reference.Channel = referenceChannelPattern) %>%
      select(Run, Stripped.Sequence.Charge, Protein.Group, Reference.Channel) -> filterSet
    
  } else {
    stop('Wrong referenceChannel parameter. Must be "H", "L", "M" or "highest"')
  }
  
  # Extract relevant columns from data based on useCols
  data %>%
    select(Run, Protein.Group, Precursor.Id, Stripped.Sequence.Charge, 
           !!sym(filterPGglobal), !!sym(filterChannel), one_of(useCols)) -> data
  
  # Summarize filtered channels based on CalcCols setting
  data %>%
    inner_join(filterSet) %>%
    filter(
      (!!sym(filterPGglobal)) < 0.01,
      (!!sym(filterChannel)) < filterChannelTh
    ) %>%
    filter_at(vars(useCols), if (CalcCols != "any") all_vars(!is.na(.) & . > 0) else any_vars(!is.na(.) & . > 0)) %>%
    group_by(Run, Stripped.Sequence.Charge, Protein.Group, Reference.Channel) %>%
    summarize(Channels.Passed.Standard.Filter = toString(str_extract(Precursor.Id, searchChannel)), .groups = "drop") -> filterSet
  
  # Flagging and summarizing for protein level
  filterSet %>%
    group_by(Run, Protein.Group) %>%
    summarize(N.Total.Precursors.Per.Protein.Standard = n(), .groups = "drop") -> protein_summary
  
  filterSet %>%
    filter(str_count(Channels.Passed.Standard.Filter, str_remove(searchChannel, "-")) >= numberChannels) %>%
    group_by(Run, Protein.Group) %>%
    summarize(N.Precursors.With.Complete.SILAC.Channels.Per.Protein.Standard = n(), .groups = "drop") -> complete_channels_summary
  
  # Combine summaries into final filterSet output
  filterSet %>%
    left_join(protein_summary, by = c("Run", "Protein.Group")) %>%
    left_join(complete_channels_summary, by = c("Run", "Protein.Group")) %>%
    select(Run, Protein.Group, Stripped.Sequence.Charge, Reference.Channel, 
           Channels.Passed.Standard.Filter, N.Total.Precursors.Per.Protein.Standard, 
           N.Precursors.With.Complete.SILAC.Channels.Per.Protein.Standard) -> filterSet
  
  # Return only relevant columns if single-channel mode
  if (numberChannels == 1) {
    filterSet %>%
      select(Run, Protein.Group, Stripped.Sequence.Charge, Reference.Channel, N.Total.Precursors.Per.Protein.Standard) -> filterSet
  }
  
  return(filterSet)
}


#' Calculation of SILAC ratios
#' 
#' Function calculating the SILAC ratios for all precursors and proteins using one of the filter sets. WARNING: Tripple SILAC CURRENTLY ONLY WORKS WITH REQUANTIFY BUT NOT FOR BASIC FILTERING. It filters for nPrecursors and then calculates all ratios anyways as long as there is any M / L (-> also if a M never passed filters)
#' 
#' @param data Output from clean_DIANN, filter_DIANN or raw DIANN report
#' @param filterSet Output from filter_DIANN
#' @param useFilter one of a) "requantify", or b) "standard". Uses Precursors for calculation where a) the reference channel passed standard filters or b) all channels passed standard filters.
#' @param precursorPerProtein How many precursors are required for one protein. Default is 1. The more precursors are required, the better quantitative performance of remaining proteins gets.
#' @param globalReference One of "L", "M", "H". Defines which channel should be used as a global reference to build normalized intensities. Should be the channel that was experimentally spiked-in.
#' @param PGcalculation Selects columns to build ratios on. Best variant should be "standard" (using Ms1.Translated and Precursor.Translated). For single cell data use "sc" (using both independent of completness). Other options are only use "Ms1.Translated" or only use "Precursor.Translated", this is not recommendet.
#' @param globalCalculation which reference proteins shall be used for the calculation of the global reference? All Heavys passing the basic filtering or all heavys having a corresponding light? 
#' @param version DIA-NN version used to analyse the raw data. Default and recommended is 1.8.1. Set to 1.9 for higher versions. Recommended DIA-NN settings for 1.8.2 beta 22 or higher are QuantUMS: legacy; turn off MBR; and --channel-spec-norm
#' 
#' @returns a tibble containing following columns: Run, Protein.Group, Strupped.Sequence.Charge, Intensity.Type (Containing "Ms1.Translated" and "Precursor.Translated"), L, M, H (if existing; contain the corresponding precursor intensities of that channel given in Intensity.Type), N.Precursors (Number of precursors defining a protein, corresponding to the respective column in filterSet that was selected in useFilter), LvsH/MvsH/LvsM (Log10 Precursor ratios of the Channels), LvsH.PG/MvsH.PG/LvsM.PG (log10 Protein.Group ratios based on the median log10 precursor ratios of both intensity types per protein & run) and Global.Log10.Reference.Intensity (Median log10 precursor intensity of the reference/spike-in channel per protein, acts as a normalization factor to calculate Abundance in the other channels. To calculate protein abundance in other channels use eg. LvsH.PG + Global.Log10.Reference.Intensity). Also contains column with information of applied filter sets and global reference used.
#' @examples 
#' # for bulk
#' data <- clean_DIANN(data) # read in the data
#' filter_frame <- filter_DIANN(data) # generate filter data frame
#' 
#' silac_ratios <- calculate_SILAC_ratios(data = data, filterSet = filter_frame, useFilter = "requantify", globalReference = "H", 
#' globalCalculation = "all", PGcalculation = "standard")
#' 
#' # for sc data
#' data <- clean_DIANN(data) # read in the data
#' filter_frame <- filter_DIANN(data,  CalcCols = "any") # generate filter data frame
#' 
#' silac_ratios <- calculate_SILAC_ratios(data = data, filterSet = filter_frame, useFilter = "requantify", globalReference = "H", 
#' globalCalculation = "all", PGcalculation = "sc")
#' 
#' # for DIA-NN 1.8.2 beta 22 or higher
#'  filter_frame <- filter_DIANN(data, CalcCols = "any", ChannelFilt = 0.01, version = "1.9")
#' silac_ratios <- calculate_SILAC_ratios(data = data, filterSet = filter_frame, useFilter = "requantify", globalReference = "H", 
#' globalCalculation = "all", PGcalculation = "standard", version = "1.9")
#' 
#' @import dplyr
#' @export
calculate_SILAC_ratios <- function(data, filterSet, useFilter = c("requantify", "standard"), 
                                   precursorPerProtein = 1, globalReference = c("L", "M", "H"), globalCalculation = c("all", "set"),
                                   PGcalculation = c("standard", "any", "Ms1.Translated", "Precursor.Translated", "sc"), version = "1.8.1"){
  
  if(useFilter == "requantify"){
    useFilter <- "N.Total.Precursors.Per.Protein.Standard"
  }else if(useFilter == "standard"){
    useFilter <- "N.Precursors.With.Complete.SILAC.Channels.Per.Protein.Standard"
    useChannels <- "Channels.Passed.Standard.Filter"
  }else{
    stop("Please define useFilter parameter.")
  }
  
  # create a precursor ID that only contains sequence and charge but not SILAC state, neccessary if using raw DIANN output
  if(!"Stripped.Sequence.Charge" %in% colnames(data)){ 
    data %>% 
      mutate(Stripped.Sequence.Charge = str_remove_all(Precursor.Id, "\\(SILAC-[KR]-[LHM]\\)")) -> data 
  } 
  
  # if wanting to use all heavy precursors for global
  if(globalCalculation == "all"){
    filterSet[,c("Run", "Protein.Group","Stripped.Sequence.Charge", useFilter)] %>% 
      rename("N.Precursors" = useFilter) %>% 
      filter(N.Precursors >= precursorPerProtein) -> filterSetGlobal
  }
  
  
  numberChannels <- max(str_count(filterSet$Channels.Passed.Standard.Filter, "[LMH]")) 
  
  if(useFilter != "N.Total.Precursors.Per.Protein.Standard"){
    filterSet[,c("Run", "Protein.Group","Stripped.Sequence.Charge", useFilter, useChannels)] %>% 
      rename("N.Precursors" = useFilter) %>% 
      rename("Passed.Precursors" = useChannels) %>% 
      filter(N.Precursors >= precursorPerProtein, 
             str_count(Passed.Precursors, "\\w") >= 2) -> filterSet
  }else{
    filterSet[,c("Run", "Protein.Group","Stripped.Sequence.Charge", useFilter)] %>% 
      rename("N.Precursors" = useFilter) %>% 
      filter(N.Precursors >= precursorPerProtein) -> filterSet
  }
  
  
  
  
  if(PGcalculation == "standard"){
    if(version == "1.8.1"){
      data %>% 
        select(Run, Protein.Group, Precursor.Id, Stripped.Sequence.Charge, Ms1.Translated, Precursor.Translated) %>% 
        mutate(Channel = str_remove(str_extract(Precursor.Id, "-[LMH]"), "-")) %>%  
        select(-Precursor.Id) %>% 
        pivot_longer(cols = contains("Translated"), names_to = "Intensity.Type", values_to = "Intensity") %>% 
        group_by(Run, Channel, Stripped.Sequence.Charge) %>% 
        filter(Intensity > 0 & !is.na(Intensity)) %>%
        # remove cases where Ms1.Translated or Precursor.Translated are missing
        filter(n() == 2) %>% 
        ungroup() %>% 
        pivot_wider(id_cols = everything(), names_from = Channel, values_from = Intensity) -> data
    }else{
      data %>% 
        select(Run, Protein.Group, Precursor.Id, Stripped.Sequence.Charge, Precursor.Normalised) %>% 
        mutate(Channel = str_remove(str_extract(Precursor.Id, "-[LMH]"), "-")) %>%  
        select(-Precursor.Id) %>% 
        pivot_longer(cols = c(Precursor.Normalised), names_to = "Intensity.Type", values_to = "Intensity") %>% 
        ungroup() %>% 
        filter(Intensity > 0 & !is.na(Intensity)) %>%
        pivot_wider(id_cols = everything(), names_from = Channel, values_from = Intensity) -> data
    }
    
  }else{
    if(version == "1.8.1"){
      data %>% 
        select(Run, Protein.Group, Precursor.Id, Stripped.Sequence.Charge, Ms1.Translated, Precursor.Translated) %>% 
        mutate(Channel = str_remove(str_extract(Precursor.Id, "-[LMH]"), "-")) %>%  
        select(-Precursor.Id) %>% 
        pivot_longer(cols = contains("Translated"), names_to = "Intensity.Type", values_to = "Intensity") %>% 
        ungroup() %>% 
        filter(Intensity > 0 & !is.na(Intensity)) %>%
        pivot_wider(id_cols = everything(), names_from = Channel, values_from = Intensity) -> data
    }else{
      data %>% 
        select(Run, Protein.Group, Precursor.Id, Stripped.Sequence.Charge, Precursor.Normalised) %>% 
        mutate(Channel = str_remove(str_extract(Precursor.Id, "-[LMH]"), "-")) %>%  
        select(-Precursor.Id) %>% 
        pivot_longer(cols = c(Precursor.Normalised), names_to = "Intensity.Type", values_to = "Intensity") %>% 
        ungroup() %>% 
        filter(Intensity > 0 & !is.na(Intensity)) %>%
        pivot_wider(id_cols = everything(), names_from = Channel, values_from = Intensity) -> data
    }
    
  }
  
  
  if(globalCalculation == "all"){
    data %>% 
      inner_join(filterSetGlobal) -> dataGlobal
  }
  
  data %>% 
    inner_join(filterSet) -> data
  
  
  if(PGcalculation %in% c("Ms1.Translated", "Precursor.Translated")){
    data %>% 
      filter(Intensity.Type == PGcalculation) -> data
    if(globalCalculation == "all"){
      dataGlobal %>% 
        filter(Intensity.Type == PGcalculation) -> dataGlobal
    }
  }
  # calculate ratios depending on the channels
  if(numberChannels == 2){
    if("H" %in% colnames(data) & "L" %in% colnames(data)){
      data %>% 
        mutate(LvsH = log10(L/H)) %>% 
        group_by(Run, Protein.Group) %>% 
        # PG level L/H ratio and L abundance normalized to global H
        mutate(LvsH.PG = median(LvsH, na.rm = T)) %>% 
        distinct() -> PG.Ratios
    }else if("H" %in% colnames(data) & "M" %in% colnames(data)){
      data %>% 
        mutate(MvsH = log10(M/H)) %>% 
        group_by(Run, Protein.Group) %>% 
        # PG level L/H ratio and L abundance normalized to global H
        mutate(MvsH.PG = median(MvsH, na.rm = T)) %>% 
        distinct() -> PG.Ratios
    }else{
      data %>% 
        mutate(LvsM = log10(L/M)) %>% 
        group_by(Run, Protein.Group) %>% 
        # PG level L/H ratio and L abundance normalized to global H
        mutate(LvsM.PG = median(LvsM, na.rm = T)) %>% 
        distinct() -> PG.Ratios
    }
    
    
    
    
    
  }else if(numberChannels == 3){
    data %>% 
      mutate(LvsH = log10(L/H), MvsH = log10(M/H), LvsM = log10(L/M)) %>% 
      group_by(Run, Protein.Group) %>% 
      # PG level L/H ratio and L abundance normalized to global H
      mutate(LvsH.PG = median(LvsH, na.rm = T), 
             MvsH.PG = median(MvsH, na.rm = T), 
             LvsM.PG = median(LvsM, na.rm = T))  %>% 
      distinct() -> PG.Ratios
  }else{
    stop("Invalid Number of Channels")
  }
  
  
  if(globalCalculation == "all"){
    dataGlobal %>% 
      ungroup %>% 
      group_by(Protein.Group, Run) %>% 
      summarize(Global.Log10.Reference.Protein.Intensity = sum(((!!sym(globalReference))), na.rm = T)) %>% 
      group_by(Protein.Group) %>% 
      #      summarize(Global.Log10.Reference.Protein.Intensity = median(log10((!!sym(globalReference))), na.rm = T)) %>% 
      summarize(Global.Log10.Reference.Protein.Intensity = median(log10(Global.Log10.Reference.Protein.Intensity), na.rm = T)) %>% 
      # mutate(Global.Log10.Reference.Protein.Intensity = median(log10((!!sym(globalReference))))) %>% 
      ungroup() %>% 
      inner_join(PG.Ratios) -> PG.Ratios
    
  }else{
    # add global protein intensity based on reference channel
    PG.Ratios %>% 
      ungroup %>% 
      group_by(Protein.Group, Run) %>% 
      summarize(Global.Log10.Reference.Protein.Intensity = sum(((!!sym(globalReference))), na.rm = T)) %>% 
      group_by(Protein.Group) %>% 
      #      summarize(Global.Log10.Reference.Protein.Intensity = median(log10((!!sym(globalReference))), na.rm = T)) %>% 
      summarize(Global.Log10.Reference.Protein.Intensity = median(log10(Global.Log10.Reference.Protein.Intensity), na.rm = T)) %>% 
      # mutate(Global.Log10.Reference.Protein.Intensity = median(log10((!!sym(globalReference))))) %>% 
      ungroup -> PG.Ratios
  }  
  
  
  
  if("NA" %in% colnames(PG.Ratios)){
    PG.Ratios %>% 
      select(-`NA`) -> PG.Ratios
  }
  
  PG.Ratios %>% 
    add_column(use.Filter = useFilter) %>% 
    add_column(precursor.Per.Protein = precursorPerProtein) %>% 
    add_column(global.Reference = globalReference) %>% 
    mutate(across(contains(".PG"), .fns=~.+Global.Log10.Reference.Protein.Intensity, .names = "{.col}.Intensity")) -> PG.Ratios
  
  
  return(PG.Ratios)
}



#' Generate human readable PG matrix outputs
#' 
#' Generates wide table versions of protein abundance and quality measurements (number of unique peptides in total as well as per run that do pass the filter). The output of this function is no input for any other function of this pipeline but to enhance easier compatibility with scripts of other users.
#' 
#' @param PGdata Output from calculate_SILAC_ratios
#' @param filterSet Output from filter_DIANN
#' 
#' @returns Returns Two wide tables: a) pg_matrix: a table with the same format as pg_matrix.tsv output from DIA-NN containing the protein abundances per run and b) peptide_counts: a wide table version of the filter_DIANN output to allow for easily readable quality control. Contains columns Protein.Group, Unique.Peptides (unique stripped peptide sequences without charge state) in all runs and L;M;H peptide counts per run. Differently from the output of filter_DIANN, this dataframe contains peptides, not precursors (no charge state).
#' @examples 
#' # for bulk
#' data <- clean_DIANN(data) # read in the data
#' filter_frame <- filter_DIANN(data) # generate filter data frame
#' 
#' silac_ratios <- calculate_SILAC_ratios(data = data, filterSet = filter_frame, useFilter = "requantify", globalReference = "H", 
#' globalCalculation = "all", PGcalculation = "standard")
#' 
#' pg_matrix <- generate_pgMatrix(silac_ratios, filter_frame)
#' 
#' @import dplyr
#' @export
generate_pgMatrix <- function(PGdata, filterSet){
  
  PGdata %>% 
    ungroup() %>% 
    select(Run, Protein.Group, contains(".PG.Intensity")) %>% 
    distinct() %>% 
    pivot_wider(id_cols = everything(), names_from = Run, values_from = contains(".PG.Intensity"), names_sep = "-", names_prefix = "Abundance.") -> PG.Intensities
  
  
  # previous counts and stuff are based on precursors, not peptides. Need to count and merge differently here
  filterSet %>% 
    mutate(Stripped.Sequence = str_remove_all(str_remove_all(Stripped.Sequence.Charge, "\\d"), "\\(UniMod:\\)")) -> filterSet
  
  filterSet %>% 
    ungroup %>% 
    filter(str_detect(Channels.Passed.Standard.Filter, "[HLM]")) %>% 
    select(Protein.Group, Stripped.Sequence) %>% 
    distinct() %>% 
    group_by(Protein.Group) %>% 
    summarise(Unique.Peptides = n()) %>% 
    full_join(filterSet) -> filterSet
  
  filterSet %>% 
    group_by(Run, Protein.Group) %>% 
    mutate(Peptide.Count.Run.L = str_count(Channels.Passed.Standard.Filter, "L"), 
           Peptide.Count.Run.M = str_count(Channels.Passed.Standard.Filter, "M"), 
           Peptide.Count.Run.H = str_count(Channels.Passed.Standard.Filter, "H")) %>% 
    select(-Stripped.Sequence.Charge, -Channels.Passed.Standard.Filter) %>%
    # removes all duplicated rows (keeps first hit)
    filter(!duplicated(Stripped.Sequence)) -> filterSet.Single
  
  # check whether in the thrown out columns has been some additional information
  filterSet %>% 
    group_by(Run, Protein.Group) %>% 
    mutate(Peptide.Count.Run.L = str_count(Channels.Passed.Standard.Filter, "L"), 
           Peptide.Count.Run.M = str_count(Channels.Passed.Standard.Filter, "M"), 
           Peptide.Count.Run.H = str_count(Channels.Passed.Standard.Filter, "H")) %>% 
    select(-Stripped.Sequence.Charge, -Channels.Passed.Standard.Filter) %>%
    # takes all rows after the first hit
    filter(duplicated(Stripped.Sequence)) -> filterSet.Doublets
  
  # if in any of the double cases (different charged peptides) a LMH was detected but not in the single cases, set 0 to 1 
  filterSet.Single %>% 
    full_join(filterSet.Doublets, by = c("Protein.Group", "Run", "Stripped.Sequence", "Unique.Peptides", 
                                         "Reference.Channel", "N.Precursors.With.Complete.SILAC.Channels.Per.Protein.Standard", 
                                         "N.Total.Precursors.Per.Protein.Standard"))  %>%
    mutate(Peptide.Count.Run.L.x = ifelse(Peptide.Count.Run.L.x == 0 & Peptide.Count.Run.L.y == 1, 1, Peptide.Count.Run.L.x), 
           Peptide.Count.Run.M.x = ifelse(Peptide.Count.Run.M.x == 0 & Peptide.Count.Run.M.y == 1, 1, Peptide.Count.Run.M.x), 
           Peptide.Count.Run.H.x = ifelse(Peptide.Count.Run.H.x == 0 & Peptide.Count.Run.H.y == 1, 1, Peptide.Count.Run.H.x)) %>% 
    ungroup %>% 
    select(Protein.Group, Run, Unique.Peptides, Stripped.Sequence, Peptide.Count.Run.L.x, Peptide.Count.Run.M.x, Peptide.Count.Run.H.x) %>% 
    distinct() %>% 
    group_by(Run, Protein.Group, Unique.Peptides) %>% 
    summarise(Peptide.Count.Run.L = sum(Peptide.Count.Run.L.x, na.rm = T), 
              Peptide.Count.Run.M = sum(Peptide.Count.Run.M.x, na.rm = T), 
              Peptide.Count.Run.H = sum(Peptide.Count.Run.H.x, na.rm = T)) %>% 
    unite(Peptide.Count, Peptide.Count.Run.L, Peptide.Count.Run.M, Peptide.Count.Run.H, sep = ";" ) %>% 
    pivot_wider(id_cols = everything(), names_from = Run, values_from = Peptide.Count, names_prefix = "Peptide.Count.") -> peptide_counts
  
  
  
  out <- list("pg_matrix" = PG.Intensities, "peptide_counts" = peptide_counts)
  return(out)
  
}


#' Sanity Check Requantify across sample ratios
#' 
#' Function to find sane requantify across sample ratios, strongly advoced to run before calculating across sample ratios
#' NOT WORKING WITH TRIPPLE SILAC
#' 
#' @param filterSet output from filter_DIANN
#' @param metadata tibble containing columns Run, Sample and replicate. Sanity Check will than be only performed across samples (replicate wise). If no metadata given, all possible Run combinations are tested
#' @param precursorPerProtein precursors per protein required, should be the same as used in calculate_SILAC_ratios
#' @returns tibble with all possible valid combinations.
#' @examples 
#' example uses bulk data
#' data <- clean_DIANN(data) # read in the data
#' filter_frame <- filter_DIANN(data) # generate filter data frame
#' 
#' silac_ratios <- calculate_SILAC_ratios(data = data, filterSet = filter_frame, useFilter = "requantify", globalReference = "H", 
#' globalCalculation = "all", PGcalculation = "standard") # calculate SILAC abundances
#' 
#' sane_PGs <- requantify_sanity_check(filterSet = filter_frame) # sanity check requantify combinations (without metadata)
#' @import dplyr
#' @export
requantify_sanity_check <- function(filterSet, metadata = NULL, precursorPerProtein = 1){
  
  # check whether a PG is valid (as in: detected in H & L channels instead of by requantify)
  
  filterSet %>%
    ungroup() %>% 
    select(Run, Protein.Group, N.Precursors.With.Complete.SILAC.Channels.Per.Protein.Standard, N.Total.Precursors.Per.Protein.Standard) %>%
    distinct() %>% 
    mutate(is.valid = ifelse(is.na(N.Total.Precursors.Per.Protein.Standard) | is.na(N.Precursors.With.Complete.SILAC.Channels.Per.Protein.Standard) |
                               N.Precursors.With.Complete.SILAC.Channels.Per.Protein.Standard < precursorPerProtein | N.Total.Precursors.Per.Protein.Standard < precursorPerProtein, 
                             0, 1)) %>% 
    select(Run, Protein.Group, is.valid) -> filterSet
  
  
  if(is.null(metadata)){
    # get matrix with all possible combinations of runs
    combinations <- combn(filterSet$Run %>% unique, 2)  
    
    # make filterset a wide tibble with Runs as column names
    filterSet %>% 
      pivot_wider(id_cols = everything(), names_from = Run, values_from = is.valid) -> filterSet
  }else if("data.frame" %in% class(metadata)){
    
    
    combinations <- combn(metadata$Sample %>% unique, 2)
    
    filterSet %>% 
      inner_join(metadata, by = "Run") %>%     
      select(-Run) %>% 
      pivot_wider(id_cols = everything(), names_from = Sample, values_from = is.valid) -> filterSet
    
    
    
    
  }else{
    warning("Object metadata is not data.frame or tibble. Run combinations will be used instead.")
    # get matrix with all possible combinations of runs
    combinations <- combn(filterSet$Run %>% unique, 2)  
    
    # make filterset a wide tibble with Runs as column names
    filterSet %>% 
      pivot_wider(id_cols = everything(), names_from = Run, values_from = is.valid) -> filterSet
  }
  
  
  
  # for a given combination do
  for( i in 1:ncol(combinations)) {
    # generate a column name that contain the combination
    filterSetColnames <- paste(combinations[1, i], combinations[2, i], sep=";")
    # genarate a new column in filterSet with the name that contains the sum of Run1 and Run2
    filterSet[, filterSetColnames] <- filterSet[, combinations[1, i]] + filterSet[, combinations[2, i]]
    # if the sum is higher than 1, the combination is valid
    
  }
  
  if(!is.null(metadata)){
    filterSet %>% 
      select(Protein.Group, replicate, contains(";")) %>%
      pivot_longer(cols = contains(";"), names_to = "Versus", values_to = "is.valid") %>% 
      filter(!is.na(is.valid) & is.valid > 0) -> validRatios
  }else{
    filterSet %>% 
      select(Protein.Group, Run, contains(";")) %>%
      pivot_longer(cols = contains(";"), names_to = "Versus", values_to = "is.valid") %>% 
      filter(!is.na(is.valid) & is.valid > 0) -> validRatios
  }
  
  
  
  return(validRatios)
  
}


#' Calculate across sample ratios
#' 
#' Function to calculate across sample ratios
#' 
#' @param data output from calculate_SILAC_ratios. In case you did something to this dataframe beforehand (eg filtering) make sure it is not grouped tibble before passing it to this function.
#' @param metadata tibble containing columns Run, Sample and replicate. Ratio calculation will than be only performed across samples (replicate wise). If no metadata given, all possible Run combinations are calculated
#' @param validRatios output from requantify_sanity_check, only neccessary if requantify was used
#' @param abundanceCol Identifier of the abundance column(s). Default is "PG.Intensity" which is the PG intensity output identifier from the calculate_SILAC_ratios function
#' @returns tibble with all across sample ratios. If validRatios was given it contains a column "is.valid". is.valid == T, the ratios is valid to use. is.valid == F ratios are returned for experimental reasons but should be removed for further analysis. 
#' @examples 
#' 
#' # example uses bulk data
#' 
#' ## for requantify data
#' data <- clean_DIANN(data) # read in the data
#' filter_frame <- filter_DIANN(data) # generate filter data frame
#' 
#' silac_ratios <- calculate_SILAC_ratios(data = data, filterSet = filter_frame, useFilter = "requantify", globalReference = "H", 
#' globalCalculation = "all", PGcalculation = "standard") # calculate SILAC abundances
#' 
#' sane_PGs <- requantify_sanity_check(filterSet = filter_frame) # sanity check requantify combinations (without metadata)
#' data_across <- calculate_across_sample_ratios(data = silac_ratios, validRatios = sanePGs)
#' 
#' 
#' ## without requantify: 
#' data <- clean_DIANN(data) # read in the data
#' filter_frame <- filter_DIANN(data) # generate filter data frame
#' 
#' silac_ratios <- calculate_SILAC_ratios(data = data, filterSet = filter_frame, useFilter = "standard", globalReference = "H", 
#' globalCalculation = "all", PGcalculation = "standard") # calculate SILAC abundances
#' data_across <- calculate_across_sample_ratios(data = silac_ratios)
#' @import dplyr
#' @export
calculate_across_sample_ratios <- function(data, metadata = NULL, validRatios = NULL, abundanceCol = "PG.Intensity"){
  
  if("list" %in% class(data)){
    data <- data[[1]]
  }
  
  if("data.frame" %in% class(metadata)){
    
    combinations <- combn(metadata$Sample %>% unique, 2)
    
    data %>% 
      select(Protein.Group, contains(abundanceCol), Run) %>% 
      distinct() %>% 
      inner_join(metadata, by = "Run") %>% 
      select(-Run) %>% 
      pivot_wider(id_cols = everything(), names_from = Sample, values_from = contains(abundanceCol)) -> data
    
  }else{
    warning("Object metadata is not data.frame or tibble. Run combinations will be used instead.")
    combinations <- combn(data$Run %>% unique, 2)  
    
    data %>% 
      select(Protein.Group, contains(abundanceCol), Run) %>% 
      distinct() %>% 
      pivot_wider(id_cols = everything(), names_from = Run, values_from = contains(abundanceCol)) -> data
  }
  
  
  # for a given combination do
  for( i in 1:ncol(combinations)) {
    # generate a column name that contain the combination
    dataColnames <- paste(combinations[1, i], combinations[2, i], sep=";")
    # genarate a new column in filterSet with the name that contains the sum of Run1 and Run2
    data[, dataColnames] <- data[, combinations[1, i]] - data[, combinations[2, i]]
    # calculate ratios of Sampole 1 and Sample 2
    
  }
  
  
  
  if(is.null(validRatios)){
    warning("No valid ratio check done, if using requantify using it is strongly adviced")
    data %>% 
      select(Protein.Group, contains(";"), contains("replicate")) %>% 
      pivot_longer(cols = contains(";"), names_to = "Versus", values_to = "Log10.Ratio") -> acrossRatios
    
  }else if("data.frame" %in% class(validRatios)){
    
    data %>% 
      select(Protein.Group, contains(";"), contains("replicate")) %>% 
      pivot_longer(cols = contains(";"), names_to = "Versus", values_to = "Log10.Ratio") %>% 
      full_join(validRatios) %>% 
      mutate(is.valid = ifelse(is.na(is.valid), F, T))-> acrossRatios
    
  }else{
    warning("Object validRatios is not data.frame or tibble. No valid ratio check done.")
    data %>% 
      select(Protein.Group, contains(";"), contains("replicate")) %>% 
      pivot_longer(cols = contains(";"), names_to = "Versus", values_to = "Log10.Ratio") -> acrossRatios
  }
  
  return(acrossRatios)
}
