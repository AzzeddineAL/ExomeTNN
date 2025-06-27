library(shiny)
library(redcapAPI)
library(tidyverse)
library(plotly)
library(patchwork)
library(DT)
library(RColorBrewer)
library(shinyjs)
library(stringdist)
library(lubridate)
library(shinydashboard)
library(edgeR)
library(limma)

options(shiny.maxRequestSize = 30 * 1024^2)

REDCap_URL <- "https://redcap.krctnn.com/api/"
REDCap_token <- "F58E14ACDF5A118B234429DC94F320F9"

SEX_LEVELS <- c("Homme", "Femme", "Non renseigné")

get_redcap_data <- function() {
  canonical_prescriber_names_list <- c(
    "laurent mesnad", "n Nicolas benichou", "nadia ouldouali", "cyril mousseaux",
    "khalil el karou", "killias bensouna", "david de saint gilles", "alexandre dupre",
    "cedric rafat", "emmanuel letavennier", "jordan lecoq", "paul gabarre",
    "hadrien debuhren", "marion delarose", "paul de lespinasse", "valentine forte",
    "armance marchal", "alexandra lehens", "jad madi", "romain brousse",
    "helene francois", "jean jacob botbol", "yannis lombardi", "justine serre",
    "ristine chenaie", "maxime chhuon", "marie imbert", "nacera ouali",
    "sylvie raballand", "rateb khayat", "caroline arches", "othmane mouhiob",
    "yves luglia", "federica bocchi", "flavie carniato", "marie lillo",
    "sarah tortonese", "xavier vincent", "coline robert", "nicolas blanchard",
    "ivan scrabine", "alexia lundy", "julie nys", "marie-camille lafargue",
    "sophie georgin lavialle", "adrien bons", "camille petit hoang", "clemence gatinolis",
    "fanny lepeytre", "paul dufour", "mathilde lemoine", "pierre galichon",
    "reda lameche", "alex lemoine", "antoine meyers", "chomroeun hov",
    "theo van den branden", "francois guy seghers", "hassan izzedine",
    "jerome tourret", "jeffrey petitguyot", "lydia arabe", "pierre antoine michel"
  )
  canonical_prescriber_names_list <- tolower(stringr::str_trim(stringr::str_squish(canonical_prescriber_names_list)))
  canonical_prescriber_names_list <- unique(canonical_prescriber_names_list)
  canonical_prescriber_names_list <- sort(canonical_prescriber_names_list)

  if (REDCap_URL == "VOTRE_URL_API_REDCAP" || REDCap_token == "VOTRE_TOKEN_API_REDCAP") {
    set.seed(123)
    data <- data.frame(
      record_id = rep(1:100, each = 2),
      ipp_exome = c(paste0("IPP", 1000 + 1:90), rep(NA, 20), paste0("IPP", 1091:1100)),
      name = paste0("Patient", rep(1:100, each = 2)),
      forname = paste0("Prenom", rep(1:100, each = 2)),
      dob_exome = sample(c(seq(as.Date('1950-01-01'), as.Date('2000-01-01'), by = "day"), NA), 200, replace = TRUE, prob = c(0.98, 0.02)),
      sex = sample(c("Homme", "Femme", "Non renseigné"), 200, replace = TRUE, prob = c(0.45, 0.45, 0.10)),
      findings_ok = sample(c("Yes", "No", "Non renseigné"), 200, replace = TRUE, prob = c(0.4, 0.5, 0.1)),
      dna_during_hospit = sample(c("Yes", "No", "Non renseigné"), 200, replace = TRUE, prob = c(0.4, 0.55, 0.05)),
      sample_date_collec = sample(c(
        seq(as.Date('2020-01-01'), as.Date('2020-12-31'), by = "day"),
        seq(as.Date('2021-01-01'), as.Date('2021-12-31'), by = "day"),
        seq(as.Date('2022-01-01'), as.Date('2022-12-31'), by = "day"),
        seq(as.Date('2023-01-01'), as.Date('2023-12-31'), by = "day"),
        seq(as.Date('2024-01-01'), Sys.Date(), by = "day"),
        Sys.Date(),
        Sys.Date() + days(1),
        NA
      ), 200, replace = TRUE, prob = c(0.1, 0.1, 0.1, 0.1, 0.3, 0.05, 0.02, 0.1)),
      emaillm = sample(c("dr.dupont@email.com", "Dr. Dupont@email.com ", "dr.martin@email.com", "Dr Martin@email.com", "DR. MARTIN@email.com", "drschmidt@email.com", "hassan.izzedine@aphp.fr", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 0.1, 0.05, 0.05)),
      first_name_referee_test = sample(c(
        "Dr. John Smith", "Jane Doe", "Prof. Alice Brown", "Charles Davis",
        "Marie Curie", "Pierre Dubois", "Sophie Martin", "Lucas Bernard",
        "HASSANE IZZEDINE", "Hassane Izzedine", "IZZEINE, H.", "HASSAN IZZED",
        "Non renseigné", NA
      ), 200, replace = TRUE, prob = c(0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.06, 0.04, 0.04, 0.04, 0.04, 0.1, 0.1)),
      origins = sample(c("Africaine", "Asiatique", "Caucasien", "Américaine", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.2, 0.15, 0.4, 0.1, 0.1, 0.05)),
      genetic_diagnosis_name = sample(c(
        "Syndrome de Rett", "Fibrose Kystique", "Drépanocytose", "Mucoviscidose",
        "Maladie de Huntington", "Syndrome de Down", "Hémophilie",
        "Maladie de Crohn", "Non renseigné", NA
      ), 200, replace = TRUE, prob = c(0.1, 0.1, 0.1, 0.1, 0.08, 0.07, 0.05, 0.05, 0.2, 0.15)),
      hgvs = sample(c("NM_000014.5:c.219G>A", "NC_000007.14:g.23345A>C", "NP_000014.2:p.Tyr23Ter", "NC_000011.10:g.54789del", "NM_000015.1:c.101T>C", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.2, 0.15, 0.1, 0.1, 0.05, 0.2, 0.2)),
      g_hgvs = sample(c("NC_000007.14:g.23345A>C", "NC_000011.10:g.54789del", "NC_000001.12:g.12345T>G", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.3, 0.2, 0.1, 0.2, 0.2)),
      p_hgvs = sample(c("NP_000014.2:p.Tyr23Ter", "NP_000015.1:p.Lys51Ile", "NP_000016.3:p.Arg123*", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.25, 0.15, 0.1, 0.3, 0.2)),
      c_hgvs = sample(c("NM_000014.5:c.219G>A", "NM_000015.1:c.101T>C", "NM_000016.3:c.367C>G", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.25, 0.15, 0.1, 0.3, 0.2)),
      variant_type = sample(c("Missense", "LoF", "Splice site", "Synonymous", "Intergenic", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.25, 0.15, 0.1, 0.1, 0.1, 0.1, 0.2)),
      gene_name = sample(c("CFTR", "BRCA1", "HTT", "SMN1", "PAH", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.2, 0.15, 0.1, 0.05, 0.05, 0.2, 0.25)),
      nanopore_test = sample(c("Yes", "No", "Non renseigné"), 200, replace = TRUE, prob = c(0.2, 0.7, 0.1)),
      wgs_compl = sample(c("Yes", "No", "Non renseigné"), 200, replace = TRUE, prob = c(0.15, 0.8, 0.05)),
      lecture1_clinic_wgslr = sample(c("OK", "À revoir", NA, "Non renseigné"), 200, replace = TRUE, prob = c(0.3, 0.2, 0.4, 0.1)),
      lecture1_clinic = sample(c("Validé", "En attente", NA, "Non renseigné"), 200, replace = TRUE, prob = c(0.3, 0.2, 0.4, 0.1)),
      rcp_conclusion = sample(c("Conclu", "Discussion nécessaire", NA, "Non renseigné"), 200, replace = TRUE, prob = c(0.3, 0.2, 0.4, 0.1)),
      rcp_completion_date = sample(c(
          seq(as.Date('2024-01-01'), Sys.Date() - days(1), by = "day"),
          Sys.Date(),
          NA
      ), 200, replace = TRUE, prob = c(0.7, 0.15, 0.15)),
      variant_rcp = sample(c("Variant A", "Variant B", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.2, 0.2, 0.5, 0.1)),
      apol1_rcp_rf = sample(c("Oui", "Non", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.1, 0.3, 0.5, 0.1)),
      rcp_inform_to_prescr = sample(c("Oui", "Non", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.4, 0.1, 0.4, 0.1)),
      rendu_results_bio = sample(c("Oui", "Non", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.3, 0.2, 0.4, 0.1)),
      incidental_findings_rcp = sample(c("Oui", "Non", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.1, 0.4, 0.4, 0.1)),
      exome_type_rcp = sample(c("0", "1", "2", "3", "Non renseigné", NA), 200, replace = TRUE, prob = c(0.2, 0.3, 0.1, 0.1, 0.2, 0.1)),
      stringsAsFactors = FALSE
    )

    data$hgvs[data$genetic_diagnosis_name == "Syndrome de Rett"][1:5] <- "NM_004992.3:c.502C>T"
    data$hgvs[data$genetic_diagnosis_name == "Fibrose Kystique"][1:5] <- "NM_000492.3:c.1521_1523delCTT"
    data$hgvs[data$genetic_diagnosis_name == "Syndrome de Rett"][1:5] <- "NM_004992.3:c.502C>T"
    data$g_hgvs[data$genetic_diagnosis_name == "Syndrome de Rett"][1:5] <- "NC_000023.11:g.153406560G>A"
    data$p_hgvs[data$genetic_diagnosis_name == "Fibrose Kystique"][1:5] <- "NP_000483.3:p.Phe508del"

    data$hgvs[data$origins == "Africaine"][1:10] <- "HGVS_Af_001"
    data$hgvs[data$origins == "Africaine"][11:20] <- "HGVS_Af_002"
    data$hgvs[data$origins == "Asiatique"][1:10] <- "HGVS_As_001"
    data$hgvs[data$origins == "Caucasien"][1:10] <- "HGVS_Cau_001"
    data$hgvs[data$origins == "Caucasien"][11:20] <- "HGVS_Cau_002"
    data$hgvs[data$origins == "Caucasien"][21:30] <- "HGVS_Cau_003"

    data$variant_type[data$hgvs == "NM_000014.5:c.219G>A"] <- "Missense"
    data$gene_name[data$hgvs == "NM_000014.5:c.219G>A"] <- "GENE_A"
    data$variant_type[data$hgvs == "NC_000007.14:g.23345A>C"] <- "LoF"
    data$gene_name[data$hgvs == "NC_000007.14:g.23345A>C"] <- "GENE_B"
    data$variant_type[data$hgvs == "NP_000014.2:p.Tyr23Ter"] <- "Splice site"
    data$gene_name[data$hgvs == "NP_000014.2:p.Tyr23Ter"] <- "GENE_C"
    data$variant_type[data$hgvs == "HGVS_Af_001"] <- "Missense"
    data$gene_name[data$hgvs == "HGVS_Af_001"] <- "GENE_AF1"
    data$variant_type[data$hgvs == "HGVS_Cau_001"] <- "LoF"
    data$gene_name[data$hgvs == "HGVS_Cau_001"] <- "GENE_CAU1"

    data_index_for_rcp_test <- which(data$record_id == 1)
    if (length(data_index_for_rcp_test) >= 2) {
      data$rcp_completion_date[data_index_for_rcp_test[1]] <- Sys.Date() - days(5)
      data$rcp_conclusion[data_index_for_rcp_test[1]] <- "Conclu"
      data$rcp_inform_to_prescr[data_index_for_rcp_test[1]] <- "Non renseigné"

      data$rcp_completion_date[data_index_for_rcp_test[2]] <- Sys.Date() + days(2)
      data$rcp_conclusion[data_index_for_rcp_test[2]] <- "Discussion nécessaire"
      data$rcp_inform_to_prescr[data_index_for_rcp_test[2]] <- "Oui"
      data$exome_type_rcp[data_index_for_rcp_test[2]] <- "3"
    }
    data_index_for_rcp_test_2 <- which(data$record_id == 2)
    if (length(data_index_for_rcp_test_2) >= 2) {
      data$rcp_completion_date[data_index_for_rcp_test_2[1]] <- Sys.Date() - days(10)
      data$rcp_conclusion[data_index_for_rcp_test_2[1]] <- "Non renseigné"

      data$rcp_completion_date[data_index_for_rcp_test_2[2]] <- Sys.Date() - days(1)
      data$rcp_conclusion[data_index_for_rcp_test_2[2]] <- "Conclu"
    }

    return(data)
  }

  tryCatch({
    rcon <- redcapAPI::redcapConnection(url = REDCap_URL, token = REDCap_token)
    data <- redcapAPI::exportRecordsTyped(rcon)

    if ("dob_exome" %in% names(data)) {
      data <- data %>%
        dplyr::mutate(age_exome = as.numeric(difftime(Sys.Date(), as.Date(dob_exome), units = "days") / 365.25))
    } else {
      data$age_exome <- NA
    }

    if (!"ipp_exome" %in% names(data)) {
      data$ipp_exome <- "Non renseigné"
    } else {
      data$ipp_exome <- as.character(data$ipp_exome)
      data$ipp_exome <- dplyr::if_else(is.na(data$ipp_exome) | data$ipp_exome == "", "Non renseigné", data$ipp_exome)
    }
    if (!"name" %in% names(data)) { data$name <- "Non renseigné" }
    if (!"forname" %in% names(data)) { data$forname <- "Non renseigné" }

    data <- data %>%
        dplyr::mutate(patient_unique_id = as.character(record_id))

    if ("sex" %in% names(data)) {
      data$sex <- as.character(data$sex)
      data$sex <- forcats::fct_explicit_na(factor(data$sex, levels = SEX_LEVELS, exclude = NULL), "Non renseigné")
    } else {
      data$sex <- factor("Non renseigné", levels = SEX_LEVELS)
    }

    if ("findings_ok" %in% names(data)) {
      data$findings_ok <- as.character(data$findings_ok)
      data$findings_ok <- forcats::fct_explicit_na(as.factor(data$findings_ok), "Non renseigné")
    } else {
      data$findings_ok <- as.factor("Non renseigné")
    }

    if ("dna_during_hospit" %in% names(data)) {
      data$dna_during_hospit <- as.character(data$dna_during_hospit)
      data$dna_during_hospit <- forcats::fct_explicit_na(as.factor(data$dna_during_hospit), "Non renseigné")
    } else {
      data$dna_during_hospit <- as.factor("Non renseigné")
    }

    if ("sample_date_collec" %in% names(data)) {
      data <- data %>%
        dplyr::mutate(sample_date_collec = as.Date(sample_date_collec)) %>%
        dplyr::mutate(collection_year = as.character(format(sample_date_collec, "%Y")),
                      collection_month = as.character(format(sample_date_collec, "%m")))
      data$collection_year <- forcats::fct_explicit_na(as.factor(data$collection_year), "Non renseigné")
      data$collection_month <- forcats::fct_explicit_na(as.factor(data$collection_month), "Non renseigné")
    } else {
      data$collection_year <- forcats::fct_explicit_na(as.factor(NA), "Non renseigné")
      data$collection_month <- forcats::fct_explicit_na(as.factor(NA), "Non renseigné")
    }

    if ("emaillm" %in% names(data)) {
      data <- data %>%
        dplyr::mutate(
          prescriber_normalized = tolower(as.character(emaillm)),
          prescriber_normalized = stringr::str_trim(prescriber_normalized),
          prescriber_normalized = stringr::str_squish(prescriber_normalized)
        )

      typo_to_canonical_map_prescriber <- c(
        "hassan.izzedine@aphp.fr" = "hassane.izzedine@aphp.fr"
      )

      original_non_renseigne_prescribers <- data$prescriber_normalized == "non renseigné" | is.na(data$prescriber_normalized) | data$prescriber_normalized == ""

      temp_prescriber_normalized <- data$prescriber_normalized
      matched_indices_prescriber <- match(temp_prescriber_normalized[!original_non_renseigne_prescribers], names(typo_to_canonical_map_prescriber))

      temp_prescriber_normalized[!original_non_renseigne_prescribers][!is.na(matched_indices_prescriber)] <- typo_to_canonical_map_prescriber[matched_indices_prescriber[!is.na(matched_indices_prescriber)]]

      data <- data %>%
        dplyr::mutate(prescriber_normalized = temp_prescriber_normalized) %>%
        dplyr::mutate(prescriber_normalized = dplyr::if_else(is.na(prescriber_normalized) | prescriber_normalized == "" | prescriber_normalized == "non renseigné", "Non renseigné", prescriber_normalized)) %>%
        dplyr::mutate(prescriber_normalized = as.factor(prescriber_normalized))
    } else {
      data$prescriber_normalized <- as.factor("Non renseigné")
    }

    if ("first_name_referee_test" %in% names(data)) {
      data <- data %>%
        dplyr::mutate(
          referee_normalized = tolower(as.character(first_name_referee_test)),
          referee_normalized = stringr::str_trim(referee_normalized),
          referee_normalized = stringr::str_squish(referee_normalized)
        ) %>%
        dplyr::mutate(referee_normalized = dplyr::if_else(
          grepl("^[0-9]+$", referee_normalized),
          "Non renseigné",
          referee_normalized
        )) %>%
        dplyr::mutate(referee_normalized = dplyr::if_else(
          is.na(referee_normalized) |
            referee_normalized == "" |
            referee_normalized == "non renseigné",
          "Non renseigné",
          referee_normalized
        ))

      typo_to_canonical_map <- c(
        "hassane izzedine" = "hassane.izzedine@aphp.fr",
        "hassan.izzedine@aphp.fr"="hassane.izzedine@aphp.fr",
        "hassane izzed" = "hassane.izzedine@aphp.fr",
        "izzedine, h." = "hassane izzedine",
        "dr izzedine"="hassane.izzedine@aphp.fr",
        "hassan izzedine"="hassane.izzedine@aphp.fr",
        "nicolas.benichou@clinique-a-pare.fr"="nicolas.benichou@aphp.fr",
        "nadia arzouk"="nadia.arzouk@aphp.fr",
        "izzedine"="hassane.izzedine@aphp.fr",
        "cyril mousseaux"="cyril.mousseaux@aphp.fr",
        "benichou"="nicolas.benichou@aphp.fr",
        "nacera ouali"="nacera.ouali@aphp.fr",
        "letavernier"="emmanuel.letavernier@aphp.fr",
        "dr. anne leblanc" = "anne-sophie leblanc",
        "anne-sophie leblanc" = "anne-sophie leblanc",
        "dr s. lebreton" = "s. lebreton",
        "dr m.l lebreton" = "s. lebreton",
        "dr pierre dubois" = "pierre dubois",
        "p.dubois" = "pierre dubois",
        "dr c.albert" = "catherine albert",
        "dr martin" = "dr. martin",
        "dr dupont" = "dr. dupont",
        "prof. alice brown" = "alice brown",
        "charles davis" = "charles davis",
        "marie curie" = "marie curie",
        "lucas bernard" = "lucas bernard",
        "john smith" = "dr. john smith",
        "jane doe" = "jane doe",
        "philippe retavergner" = "philippe retavergnier",
        "arnaud renau" = "arnauld renauld",
        "jean louis poignet" = "jean-louis poignet",
        "dr. john smith" = "dr. john smith",
        "khalil.el-karou@aphp.fr" = "khalil el karou"
      )

      non_na_and_non_numeric_refs <- data$referee_normalized[data$referee_normalized != "Non renseigné"]
      matched_indices <- match(non_na_and_non_numeric_refs, names(typo_to_canonical_map))
      data$referee_normalized[data$referee_normalized != "Non renseigné"][!is.na(matched_indices)] <- typo_to_canonical_map[matched_indices[!is.na(matched_indices)]]

      unique_referees_to_process <- unique(data$referee_normalized[data$referee_normalized != "Non renseigné"])

      referee_fuzzy_map <- setNames(unique_referees_to_process, unique_referees_to_process)

      if (length(unique_referees_to_process) > 0 && length(canonical_prescriber_names_list) > 0) {
        for (ref_name in unique_referees_to_process) {
          if (ref_name == "non renseigné") next

          distances <- stringdist::stringdist(ref_name, canonical_prescriber_names_list, method = "lv")
          min_dist <- min(distances)
          threshold <- 3

          if (min_dist <= threshold) {
            best_match_canonical <- canonical_prescriber_names_list[which.min(distances)]
            if (ref_name != best_match_canonical) {
              referee_fuzzy_map[ref_name] <- best_match_canonical
            }
          }
        }
      }
      data$referee_normalized <- dplyr::recode(data$referee_normalized, !!!referee_fuzzy_map)
      data$first_name_referee_test <- forcats::fct_explicit_na(factor(data$referee_normalized, levels = unique(data$referee_normalized), exclude = NULL), "Non renseigné")
      data <- data %>% dplyr::select(-referee_normalized)
    } else {
      data$first_name_referee_test <- as.factor("Non renseigné")
    }

    if ("origins" %in% names(data)) {
      data$origins <- as.character(data$origins)
      data$origins <- forcats::fct_explicit_na(as.factor(data$origins), "Non renseigné")
    } else {
      data$origins <- as.factor("Non renseigné")
    }

    if ("genetic_diagnosis_name" %in% names(data)) {
      data$genetic_diagnosis_name <- as.character(data$genetic_diagnosis_name)
      data$genetic_diagnosis_name <- forcats::fct_explicit_na(as.factor(data$genetic_diagnosis_name), "Non renseigné")
    } else {
      data$genetic_diagnosis_name <- as.factor("Non renseigné")
    }

    variant_cols <- c("hgvs", "g_hgvs", "p_hgvs", "c_hgvs")
    for (col in variant_cols) {
      if (col %in% names(data)) {
        data[[col]] <- as.character(data[[col]])
        data[[col]] <- dplyr::if_else(is.na(data[[col]]) | stringr::str_trim(tolower(data[[col]])) %in% c("", "non renseigné"), "Non renseigné", data[[col]])
      } else {
        data[[col]] <- "Non renseigné"
      }
    }

    if ("variant_type" %in% names(data)) {
      data$variant_type <- as.character(data$variant_type)
      data$variant_type <- forcats::fct_explicit_na(as.factor(data$variant_type), "Non renseigné")
    } else {
      data$variant_type <- as.factor("Non renseigné")
    }

    if ("gene_name" %in% names(data)) {
      data$gene_name <- as.character(data$gene_name)
      data$gene_name <- dplyr::if_else(is.na(data$gene_name) | stringr::str_trim(tolower(data$gene_name)) %in% c("", "non renseigné"), "Non renseigné", data[[col]])
    } else {
      data$gene_name <- "Non renseigné"
    }

    if ("nanopore_test" %in% names(data)) {
      data$nanopore_test <- as.character(data$nanopore_test)
      data$nanopore_test <- forcats::fct_explicit_na(as.factor(data$nanopore_test), "Non renseigné")
    } else {
      data$nanopore_test <- as.factor("Non renseigné")
    }

    if ("wgs_compl" %in% names(data)) {
      data$wgs_compl <- as.character(data$wgs_compl)
      data$wgs_compl <- forcats::fct_explicit_na(as.factor(data$wgs_compl), "Non renseigné")
    } else {
      data$wgs_compl <- as.factor("Non renseigné")
    }

    rcp_vars_to_check <- c("lecture1_clinic_wgslr", "lecture1_clinic", "rcp_conclusion",
                           "variant_rcp", "apol1_rcp_rf", "rcp_inform_to_prescr",
                           "rendu_results_bio", "incidental_findings_rcp", "exome_type_rcp")
    for (col_name in rcp_vars_to_check) {
      if (!col_name %in% names(data)) {
        data[[col_name]] <- "Non renseigné"
      } else {
        data[[col_name]] <- as.character(data[[col_name]])
        data[[col_name]] <- dplyr::if_else(is.na(data[[col_name]]) | stringr::str_trim(tolower(data[[col_name]])) %in% c("", "non renseigné"), "Non renseigné", data[[col_name]])
      }
    }
    if (!"rcp_completion_date" %in% names(data)) {
        data$rcp_completion_date <- as.Date(data$sample_date_collec)
    } else {
        data$rcp_completion_date <- suppressWarnings(as.Date(data$rcp_completion_date))
    }

    data <- data %>%
      dplyr::arrange(
        patient_unique_id,
        dplyr::desc(rcp_completion_date),
        dplyr::desc(sample_date_collec)
      ) %>%
      dplyr::group_by(patient_unique_id) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()

    return(data)
  }, error = function(e) {
    showNotification("Erreur de connexion REDCap. Vérifiez URL/Token et permissions. Utilisation de données fictives par défaut.", type = "error", duration = 8)
    data.frame(
      record_id = character(), ipp_exome = character(), name = character(), forname = character(),
      dob_exome = as.Date(character()), sex = factor(levels = SEX_LEVELS), findings_ok = factor(levels = c("Yes", "No", "Non renseigné")),
      dna_during_hospit = factor(), sample_date_collec = as.Date(character()),
      collection_year = factor(), collection_month = factor(),
      emaillm = character(), first_name_referee_test = factor(),
      origins = factor(), genetic_diagnosis_name = factor(), hgvs = character(),
      g_hgvs = character(), p_hgvs = character(), c_hgvs = character(),
      variant_type = factor(), gene_name = character(),
      nanopore_test = factor(), wgs_compl = factor(),
      lecture1_clinic_wgslr = character(),
      lecture1_clinic = character(),
      rcp_conclusion = character(),
      rcp_completion_date = as.Date(character()),
      variant_rcp = character(),
      apol1_rcp_rf = character(),
      rcp_inform_to_prescr = character(),
      rendu_results_bio = character(),
      incidental_findings_rcp = character(),
      exome_type_rcp = character(),
      patient_unique_id = character(),
      stringsAsFactors = FALSE
    )
  })
}

ui <- dashboardPage(
  dashboardHeader(title = "Dashboard REDCap Exome"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Vue d'ensemble", tabName = "overview", icon = icon("dashboard")),
      menuItem("Nouveaux Patients", tabName = "new_patients", icon = icon("user-plus")),
      menuItem("Prélèvements", tabName = "prelevements", icon = icon("syringe")),
      menuItem("Médecins Responsables", tabName = "prescribers", icon = icon("user-md")),
      menuItem("Médecins Référents", tabName = "referrers", icon = icon("user-nurse")),
      menuItem("Diagnostics Génétiques", tabName = "diagnostics", icon = icon("dna"),
        uiOutput("genetic_diagnosis_filter_ui")
      ),

      menuItem("Variants", tabName = "variants", icon = icon("microscope"),
        menuSubItem("Variants Récurrents", tabName = "variants_recurrents", icon = icon("chart-bar")),
        menuSubItem("Variants par Maladie", tabName = "variants_by_disease", icon = icon("disease")),
        menuSubItem("Variants par Ethnie", tabName = "variants_by_ethnie", icon = icon("globe-americas"))
      ),
      menuItem("Analyse RNA-Seq", tabName = "rnaseq_analysis", icon = icon("analytics")),
      menuItem("Données brutes", tabName = "raw_data", icon = icon("table")),

      selectInput("filter_findings_ok", "Filtrer par Résultat Exome :",
        choices = c("Tous", "Tous (Excl. NA)", "Yes", "No", "Non renseigné"), selected = "Tous"
      ),
      selectInput("filter_hospit", "Filtrer par Prélèvement en Hospi :",
        choices = c("Tous", "Tous (Excl. NA)", "Yes", "No", "Non renseigné"), selected = "Tous"
      ),
      uiOutput("year_filter_ui")
    ),
    tags$div(
      style = "position: absolute; bottom: 10px; width: 100%; text-align: center; padding: 0 10px; box-sizing: border-box; color: #b8c7ce;",
      tags$p(style = "margin-bottom: 0;", "Powered with Shiny-R"),
      tags$p(style = "margin-top = 0;", "Developed by Azzeddine A.")
    )
  ),
  dashboardBody(
    shinyjs::useShinyjs(),
    tabItems(
      tabItem(
        tabName = "overview",
        h2("Statistiques Clés du Centre"),
        fluidRow(
          valueBoxOutput("total_patients_box"),
          valueBoxOutput("avg_age_box"),
          valueBoxOutput("positive_results_perc_box"),
          valueBoxOutput("total_prescribers_box"),
          valueBoxOutput("patients_without_ipp_box")
        ),
        fluidRow(
          valueBoxOutput("nanopore_seq_box"),
          valueBoxOutput("wgs_seq_box")
        ),
        fluidRow(
          box(
            title = "Répartition par Sexe", status = "primary", solidHeader = TRUE,
            collapsible = TRUE,
            plotlyOutput("sex_distribution_plot")
          ),
          box(
            title = "Distribution des Résultats d'Exomes", status = "primary", solidHeader = TRUE,
            collapsible = TRUE,
            plotlyOutput("findings_ok_distribution_plot")
          )
        ),
        fluidRow(
          box(
            title = "Distribution de l'Âge", status = "primary", solidHeader = TRUE,
            collapsible = TRUE,
            plotlyOutput("age_distribution_plot")
          ),
          box(
            title = "Distribution des Origines", status = "primary", solidHeader = TRUE,
            collapsible = TRUE,
            plotlyOutput("origins_distribution_plot")
          )
        ),
        fluidRow(
          box(
            title = "Évolution du Nombre de Patients par Année", status = "primary", solidHeader = TRUE,
            collapsible = TRUE, width = 12,
            plotlyOutput("patients_by_year_plot")
          )
        )
      ),

      tabItem(
        tabName = "new_patients",
        h2("Nouveaux Patients"),
        fluidRow(
          valueBoxOutput("new_patients_previous_month_box"),
          valueBoxOutput("new_patients_current_month_box")
        ),
        fluidRow(
          box(
            title = textOutput("patients_previous_month_title"),
            status = "success", solidHeader = TRUE,
            collapsible = TRUE, width = 6,
            DT::dataTableOutput("patients_previous_month_table"),
            downloadButton("download_patients_previous_month", "Télécharger ce tableau (CSV)")
          ),
          box(
            title = textOutput("patients_current_month_title"),
            status = "info", solidHeader = TRUE,
            collapsible = TRUE, width = 6,
            DT::dataTableOutput("patients_current_month_table"),
            downloadButton("download_patients_current_month", "Télécharger ce tableau (CSV)")
          )
        ),
        
        hr(),
        h2("Télécharger les Patients par Intervalle de Dates"),
        fluidRow(
          box(
            title = "Sélectionner une période et télécharger", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 12,
            column(6,
                   dateRangeInput(
                     "download_patients_date_range",
                     "Sélectionner l'intervalle de dates pour le prélèvement :",
                     start = Sys.Date() - years(1),
                     end = Sys.Date(),
                     format = "yyyy-mm-dd",
                     separator = "à"
                   )
            ),
            column(6,
                   br(),
                   downloadButton("download_selected_date_range_patients", "Télécharger les Patients de la Période (CSV)")
            )
          )
        ),
        
        hr(),
        h2("Suivi des Mises à Jour RCP"),
        fluidRow(
            box(
                title = "Contrôles", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 4,
                p("Les données des mises à jour RCP sont chargées automatiquement. Elles sont rafraîchies lorsque les données brutes principales sont rafraîchies."),
                
                hr(),
                
                h4("Filtrage du Journal (par date)"),
                dateRangeInput(
                    "log_rcp_date_range",
                    "Sélectionner un intervalle de dates:",
                    start = Sys.Date() - days(7),
                    end = Sys.Date(),
                    format = "yyyy-mm-dd",
                    separator = "à"
                ),
                p("Cette section affiche la mise à jour RCP la plus récente pour chaque patient dont le dossier a été modifié sur les variables RCP spécifiques."),
                uiOutput("log_rcp_status_message")
            ),
            
            box(
                title = "Mises à Jour RCP des Patients (Dernière par Patient)", status = "success", solidHeader = TRUE, collapsible = TRUE, width = 8,
                DTOutput("log_rcp_log_data_table"),
                downloadButton("log_rcp_download_data", "Télécharger le journal affiché (.csv)")
            )
        )
      ),

      tabItem(
        tabName = "prelevements",
        h2("Statistiques sur les Prélèvements (DNA)"),
        fluidRow(
          valueBoxOutput("hospit_yes_box"),
          valueBoxOutput("hospit_no_box")
        ),
        fluidRow(
          box(
            title = "Répartition des Prélèvements", status = "info", solidHeader = TRUE,
            collapsible = TRUE,
            plotlyOutput("dna_hospit_distribution_plot")
          )
        )
      ),

      tabItem(
        tabName = "prescribers",
        h2("Statistiques par Prescripteur"),
        fluidRow(
          box(
            title = "Nombre de Prescriptions par Médecin", status = "success", solidHeader = TRUE,
            collapsible = TRUE, width = 6,
            plotlyOutput("prescriptions_by_prescriber_plot")
          ),
          box(
            title = "Statistiques des Résultats Exome par Médecin", status = "success", solidHeader = TRUE,
            collapsible = TRUE, width = 6,
            DT::dataTableOutput("findings_by_prescriber_table")
          )
        ),
        fluidRow(
          box(
            title = "Détail des Prescriptions par Médecin (Top 10)", status = "success", solidHeader = TRUE,
            collapsible = TRUE, width = 12,
            plotlyOutput("top_prescribers_findings_plot")
          )
        ),
        fluidRow(
          box(
            title = textOutput("top_prescribers_current_month_title"), status = "primary", solidHeader = TRUE,
            collapsible = TRUE, width = 6,
            plotlyOutput("top_prescribers_current_month_plot")
          ),
          box(
            title = textOutput("top_prescribers_previous_month_title"), status = "primary", solidHeader = TRUE,
            collapsible = TRUE, width = 6,
            plotlyOutput("top_prescribers_previous_month_plot")
          )
        )
      ),

      tabItem(
        tabName = "referrers",
        h2("Statistiques par Médecin Référent"),
        fluidRow(
          box(
            title = "Nombre de Patients Adressés par Médecin Référent", status = "primary", solidHeader = TRUE,
            collapsible = TRUE, width = 6,
            plotlyOutput("referrals_by_referrer_plot")
          ),
          box(
            title = "Statistiques des Résultats Exome par Médecin Référent", status = "primary", solidHeader = TRUE,
            collapsible = TRUE, width = 6,
            DT::dataTableOutput("findings_by_referrer_table")
          )
        ),
        fluidRow(
          box(
            title = "Détail des Résultats Exome par Référent (Top 10)", status = "primary", solidHeader = TRUE,
            collapsible = TRUE, width = 12,
            plotlyOutput("top_referrers_findings_plot")
          )
        ),
        fluidRow(
          box(
            title = textOutput("top_referrers_current_month_title"), status = "info", solidHeader = TRUE,
            collapsible = TRUE, width = 6,
            plotlyOutput("top_referrers_current_month_plot")
          ),
          box(
            title = textOutput("top_referrers_previous_month_title"), status = "info", solidHeader = TRUE,
            collapsible = TRUE, width = 6,
            plotlyOutput("top_referrers_previous_month_plot")
          )
        )
      ),

      tabItem(
        tabName = "diagnostics",
        h2("Patients par Diagnostic Génétique"),
        fluidRow(
          box(
            title = "Patients Atteints du Diagnostic Sélectionné", status = "info", solidHeader = TRUE,
            collapsible = TRUE, width = 12,
            DT::dataTableOutput("patients_by_diagnosis_table")
          )
        )
      ),

      tabItem(
        tabName = "variants_recurrents",
        h2("Analyse des Variants Génétiques - Variants Récurrents (Global)"),
        fluidRow(
          valueBoxOutput("total_variants_analyzed_box"),
          valueBoxOutput("hgvs_na_variants_box"),
          valueBoxOutput("top_hgvs_variant_box")
        ),
        box(
          title = "Variants HGVS les Plus Fréquents (>1% des variants renseignés)", status = "primary", solidHeader = TRUE,
          collapsible = TRUE, width = 12,
          plotlyOutput("recurring_variants_plot")
        )
      ),
      tabItem(
        tabName = "variants_by_disease",
        h2("Analyse des Variants Génétiques - Variants par Maladie"),
        fluidRow(
          box(
            title = "Sélectionner la/les Maladie(s)", status = "info", solidHeader = TRUE,
            width = 4,
            uiOutput("disease_filter_ui_variants")
          ),
          box(
            title = "Détails des Variants pour la Maladie Sélectionnée", status = "info", solidHeader = TRUE,
            width = 8,
            DT::dataTableOutput("variants_by_disease_table")
          )
        )
      ),
      tabItem(
        tabName = "variants_by_ethnie",
        h2("Analyse des Variants Génétiques - Variants par Ethnie"),
        fluidRow(
          box(
            title = "Sélectionner l'Ethnie (ou les Ethnies)", status = "warning", solidHeader = TRUE,
            width = 4,
            uiOutput("origins_filter_ui_variants")
          )
        ),
        fluidRow(
          box(
            title = "Top Variant Protéique (p_hgvs) par Ethnie", status = "info", solidHeader = TRUE,
            collapsible = TRUE, width = 4,
            plotlyOutput("top_p_hgvs_by_ethnie_plot")
          ),
          box(
            title = "Top Variant Codant (c_hgvs) par Ethnie", status = "info", solidHeader = TRUE,
            collapsible = TRUE, width = 4,
            plotlyOutput("top_c_hgvs_by_ethnie_plot")
          ),
          box(
            title = "Top Variant Génomique (g_hgvs) par Ethnie", status = "info", solidHeader = TRUE,
            collapsible = TRUE, width = 4,
            plotlyOutput("top_g_hgvs_by_ethnie_plot")
          )
        )
      ),

      tabItem(
        tabName = "rnaseq_analysis",
        h2("Analyse d'Expression Différentielle RNA-Seq (edgeR)"),
        fluidRow(
          box(
            title = "Chargement des Données et Paramètres d'Analyse", status = "primary", solidHeader = TRUE, width = 12,
            column(
              6,
              fileInput("counts_file", "Choisir la matrice de comptage (CSV/TSV)",
                multiple = FALSE,
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  "text/tab-separated-values",
                  ".tsv",
                  ".txt"
                )
              ),
              downloadButton("download_example_counts", "Télécharger exemple matrice de comptage (CSV)")
            ),
            column(
              6,
              fileInput("metadata_file", "Choisir le fichier de métadonnées (CSV/TSV)",
                multiple = FALSE,
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv",
                  "text/tab-separated-values",
                  ".tsv",
                  ".txt"
                )
              ),
              downloadButton("download_example_metadata", "Télécharger exemple métadonnées (CSV)")
            ),
            hr(),
            uiOutput("condition_selector_ui"),
            uiOutput("batch_selector_ui"),
            numericInput("logfc_threshold", "Seuil de Log2 Fold Change (logFC) :", value = 1, min = 0, step = 0.1),
            numericInput("fdr_threshold", "Seuil de FDR (False Discovery Rate) :", value = 0.05, min = 0, max = 1, step = 0.01),
            actionButton("run_rnaseq_analysis", "Lancer l'Analyse RNA-Seq", icon = icon("play")),
            br(), br(),
            verbatimTextOutput("rnaseq_status")
          )
        ),
        fluidRow(
          box(
            title = "Visualisations de l'Analyse RNA-Seq", status = "info", solidHeader = TRUE, width = 12,
            tabsetPanel(
              tabPanel("MDS Plot", plotlyOutput("mds_plot")),
              tabPanel("MA Plot", plotlyOutput("ma_plot")),
              tabPanel("Volcano Plot", plotlyOutput("volcano_plot")),
              tabPanel("Heatmap des Distances Échantillon-Échantillon", plotlyOutput("heatmap_plot"))
            )
          )
        ),
        fluidRow(
          box(
            title = "Résultats de l'Analyse RNA-Seq (Gènes Différentiellement Exprimés)", status = "success", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("rnaseq_results_table"),
            downloadButton("download_rnaseq_results", "Télécharger les Résultats (CSV)")
          )
        )
      ),

      tabItem(
        tabName = "raw_data",
        h2("Données Brutes (Premières Lignes)"),
        actionButton("refresh_raw_data", "Rafraîchir les données brutes"),
        hr(),
        DT::dataTableOutput("raw_data_table")
      )
    )
  )
)

server <- function(input, output, session) {

  options(shiny.fullstacktrace = TRUE)

  rnaseq_analysis_data <- reactiveValues(
    counts = NULL,
    metadata = NULL,
    dge = NULL,
    lrt = NULL,
    results_table = NULL,
    status = "En attente du chargement des fichiers..."
  )

  read_file_data <- function(filepath, file_ext, type = c("counts", "metadata")) {
    type <- match.arg(type)
    df_raw <- NULL

    seps_to_try <- character(0)
    if (file_ext == "csv") {
      seps_to_try <- c(",", ";", "\t")
    } else if (file_ext == "tsv" || file_ext == "txt") {
      seps_to_try <- c("\t", ",", ";")
    } else {
      seps_to_try <- c("\t", ",", ";")
    }

    found_multi_column_df <- FALSE
    for (s in seps_to_try) {
      tryCatch({
        temp_df <- suppressWarnings(utils::read.delim(filepath, sep = s, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, fill = TRUE))

        if (ncol(temp_df) > 1) {
          df_raw <- temp_df
          found_multi_column_df <- TRUE
          break
        }
      }, error = function(e) {
      })
    }

    if (is.null(df_raw) || nrow(df_raw) == 0) {
      stop("Le fichier est vide ou n'a pas pu être lu avec les séparateurs courants (virgule, point-virgule, tabulation). Assurez-vous qu'il contient des données et des en-têtes valides.")
    }

    if (!found_multi_column_df && ncol(df_raw) == 1) {
        first_col_name_final <- names(df_raw)[1]
        if (type == "counts") {
            stop(paste("Le fichier de comptage ne contient qu'une seule colonne d'identifiants et aucune donnée de comptage. Le fichier doit contenir au moins une colonne de comptage. Vérifiez le séparateur (virgule, point-virgule, tabulation) ou l'en-tête, et assurez-vous que le fichier n'est pas une simple liste d'identifiants.", sep=""))
        } else if (type == "metadata") {
            stop(paste("Le fichier de métadonnées ne contient qu'une seule colonne ('", first_col_name_final, "') qui semble être des identifiants et aucune information de condition. Le fichier doit contenir au moins une colonne supplémentaire pour les conditions expérimentales. Vérifiez le séparateur (virgule, point-virgule, tabulation) ou l'en-tête.", sep=""))
        }
    }

    first_col_name <- names(df_raw)[1]
    if (first_col_name == "") {
        names(df_raw)[1] <- "ID_col"
        first_col_name <- "ID_col"
    }

    ids <- df_raw[[first_col_name]]

    if (any(duplicated(ids))) {
        stop(paste("Identifiants dupliqués détectés dans la première colonne ('", first_col_name, "'). Veuillez vous assurer que tous les identifiants de gènes/échantillons sont uniques.", sep=""))
    }

    df <- df_raw[, -which(names(df_raw) == first_col_name), drop = FALSE]

    rownames(df) <- ids

    if (type == "counts") {
        for (col_idx in 1:ncol(df)) {
            col_name_current <- names(df)[col_idx]
            original_col <- df[[col_name_current]]
            numeric_col <- suppressWarnings(as.numeric(as.character(original_col)))

            non_numeric_elements_idx <- which(is.na(numeric_col) & !is.na(original_col))

            if (length(non_numeric_elements_idx) > 0) {
                example_non_numeric_value <- original_col[non_numeric_elements_idx[1]]
                stop(paste0("La matrice de comptage doit contenir uniquement des valeurs numériques. Des valeurs non-numériques ont été détectées dans la colonne '", col_name_current, "' (ex: '", example_non_numeric_value, "'). Assurez-vous qu'il n'y a pas de texte ou de caractères spéciaux dans les comptages, et que la première colonne est bien les identifiants de gènes/transcrits."))
            }
        }
        df <- as.data.frame(lapply(df, as.numeric))

        if (any(df < 0) || any(df %% 1 != 0)) {
            stop("Les comptages doivent être des entiers non-négatifs.")
        }
    } else if (type == "metadata") {
        names(df) <- make.names(names(df), unique = TRUE)
    }
    return(df)
  }

  observeEvent(input$counts_file, {
    req(input$counts_file)
    tryCatch(
      {
        file_ext <- tools::file_ext(input$counts_file$datapath)
        rnaseq_analysis_data$counts <- read_file_data(input$counts_file$datapath, file_ext, "counts")
        rnaseq_analysis_data$status <- "Matrice de comptage chargée. Chargez les métadonnées."
      },
      error = function(e) {
        rnaseq_analysis_data$counts <- NULL
        error_msg <- if (!is.null(e$message)) as.character(e$message) else "Erreur inconnue."
        rnaseq_analysis_data$status <- paste("Erreur de chargement de la matrice de comptage:", error_msg)
        showNotification("Erreur de chargement de la matrice de comptage.", type = "error")
      }
    )
  })

  observeEvent(input$metadata_file, {
    req(input$metadata_file)
    tryCatch(
      {
        file_ext <- tools::file_ext(input$metadata_file$datapath)
        rnaseq_analysis_data$metadata <- read_file_data(input$metadata_file$datapath, file_ext, "metadata")
        rnaseq_analysis_data$status <- "Fichier de métadonnées chargé. Vous pouvez lancer l'analyse."
      },
      error = function(e) {
        rnaseq_analysis_data$metadata <- NULL
        error_msg <- if (!is.null(e$message)) as.character(e$message) else "Erreur inconnue."
        rnaseq_analysis_data$status <- paste("Erreur de chargement du fichier de métadonnées:", error_msg)
        showNotification("Erreur de chargement du fichier de métadonnées.", type = "error")
      }
    )
  })

  output$download_example_counts <- downloadHandler(
    filename = function() {
      "exemple_matrice_comptage.csv"
    },
    content = function(file) {
      example_counts <- data.frame(
        GeneID = paste0("Gene", 1:10),
        SampleA_Rep1 = sample(100:1000, 10),
        SampleA_Rep2 = sample(100:1000, 10),
        SampleB_Rep1 = sample(100:1000, 10),
        SampleB_Rep2 = sample(100:1000, 10)
      )
      write.csv(example_counts, file, row.names = FALSE)
    }
  )

  output$download_example_metadata <- downloadHandler(
    filename = function() {
      "exemple_metadonnees.csv"
    },
    content = function(file) {
      example_metadata <- data.frame(
        SampleID = c("SampleA_Rep1", "SampleA_Rep2", "SampleB_Rep1", "SampleB_Rep2"),
        Condition = c("Control", "Control", "Treated", "Treated"),
        Batch = c("Batch1", "Batch2", "Batch1", "Batch2")
      )
      write.csv(example_metadata, file, row.names = FALSE)
    }
  )

  output$condition_selector_ui <- renderUI({
    req(rnaseq_analysis_data$metadata)
    valid_cols <- names(rnaseq_analysis_data$metadata)[sapply(rnaseq_analysis_data$metadata, function(x) !is.numeric(x) || length(unique(x)) < nrow(rnaseq_analysis_data$metadata) * 0.5)]
    if (length(valid_cols) == 0) {
      return(tagList(
        p("Aucune colonne appropriée pour la condition trouvée dans les métadonnées."),
        p("Veuillez vous assurer que vos colonnes ne sont pas toutes numériques ou ont trop de valeurs uniques.")
      ))
    }
    selectInput("rnaseq_condition_col", "Choisir la colonne de Condition :", choices = valid_cols)
  })

  output$batch_selector_ui <- renderUI({
    req(rnaseq_analysis_data$metadata)
    valid_cols <- names(rnaseq_analysis_data$metadata)[sapply(rnaseq_analysis_data$metadata, function(x) !is.numeric(x) || length(unique(x)) < nrow(rnaseq_analysis_data$metadata) * 0.5)]
    selectInput("rnaseq_batch_col", "Choisir la colonne de Lot (Batch - optionnel) :",
      choices = c("Aucun", valid_cols), selected = "Aucun"
    )
  })

  output$rnaseq_status <- renderText({
    rnaseq_analysis_data$status
  })

  observeEvent(input$run_rnaseq_analysis, {
    req(rnaseq_analysis_data$counts, rnaseq_analysis_data$metadata, input$rnaseq_condition_col)

    rnaseq_analysis_data$status <- "Validation des données..."

    counts <- rnaseq_analysis_data$counts
    metadata <- rnaseq_analysis_data$metadata
    condition_col <- input$rnaseq_condition_col
    batch_col <- input$rnaseq_batch_col

    common_samples <- intersect(colnames(counts), rownames(metadata))
    if (length(common_samples) == 0) {
      rnaseq_analysis_data$status <- "Erreur: Aucun nom d'échantillon commun entre la matrice de comptage et les métadonnées."
      showNotification(rnaseq_analysis_data$status, type = "error", duration = 8)
      return()
    }
    if (length(common_samples) < ncol(counts) || length(common_samples) < nrow(metadata)) {
      showNotification("Attention: Certains échantillons ne sont pas présents dans les deux fichiers et seront ignorés.", type = "warning", duration = 5)
    }

    counts <- counts[, common_samples, drop = FALSE]
    metadata <- metadata[common_samples, , drop = FALSE]

    if (!condition_col %in% names(metadata)) {
      rnaseq_analysis_data$status <- paste("Erreur: La colonne de condition '", condition_col, "' n'existe pas dans les métadonnées.", sep = "")
      showNotification(rnaseq_analysis_data$status, type = "error", duration = 8)
      return()
    }
    metadata[[condition_col]] <- factor(as.character(metadata[[condition_col]]))
    if (length(unique(metadata[[condition_col]])) < 2) {
      rnaseq_analysis_data$status <- "Erreur: La colonne de condition doit avoir au moins deux niveaux uniques pour l'analyse différentielle."
      showNotification(rnaseq_analysis_data$status, type = "error", duration = 8)
      return()
    }

    zero_count_samples <- colnames(counts)[colSums(counts) == 0]
    if (length(zero_count_samples) > 0) {
      rnaseq_analysis_data$status <- paste("Attention: Les échantillons suivants ont 0 comptages totaux et seront supprimés:", paste(zero_count_samples, collapse = ", "))
      showNotification(rnaseq_analysis_data$status, type = "warning", duration = 8)

      counts <- counts[, !colnames(counts) %in% zero_count_samples, drop = FALSE]
      metadata <- metadata[!rownames(metadata) %in% zero_count_samples, , drop = FALSE]

      if (ncol(counts) < 2 || nrow(metadata) < 2) {
        rnaseq_analysis_data$status <- "Erreur: Après suppression des échantillons à 0 comptages, il n'y a pas assez d'échantillons pour l'analyse RNA-Seq."
        showNotification(rnaseq_analysis_data$status, type = "error", duration = 8)
        return()
      }

      if (length(unique(metadata[[condition_col]])) < 2) {
        rnaseq_analysis_data$status <- "Erreur: Après suppression des échantillons à 0 comptages, la colonne de condition n'a plus assez de niveaux uniques pour l'analyse différentielle."
        showNotification(rnaseq_analysis_data$status, type = "error", duration = 8)
        return()
      }
    }

    design_formula_str <- paste("~", condition_col)
    if (batch_col != "Aucun" && batch_col %in% names(metadata)) {
      metadata[[batch_col]] <- factor(as.character(metadata[[batch_col]]))
      design_formula_str <- paste("~", batch_col, "+", condition_col)

      if (any(table(metadata[[condition_col]], metadata[[batch_col]]) == 0)) {
        showNotification("Attention: Certaines combinaisons Condition/Batch n'ont pas d'échantillons, cela peut causer des problèmes dans le design.", type = "warning", duration = 8)
      }
    }

    rnaseq_analysis_data$status <- paste("Formule de design utilisée:", design_formula_str)

    tryCatch(
      {
        design_formula <- as.formula(design_formula_str)

        if (!all(all.vars(design_formula) %in% names(metadata))) {
          stop("Certaines variables dans la formule de design ne sont pas présentes dans le fichier de métadonnées.")
        }

        rnaseq_analysis_data$status <- "Exécution de l'analyse edgeR..."
        dge <- edgeR::DGEList(counts = counts, group = metadata[[condition_col]])
        keep <- edgeR::filterByExpr(dge, design = model.matrix(design_formula, data = metadata))
        dge <- dge[keep, , keep.lib.sizes = FALSE]
        dge <- edgeR::calcNormFactors(dge)
        design_matrix <- model.matrix(design_formula, data = metadata)
        dge <- edgeR::estimateDisp(dge, design_matrix, robust = TRUE)
        fit <- edgeR::glmFit(dge, design_matrix)
        lrt <- edgeR::glmLRT(fit, coef = ncol(design_matrix))

        rnaseq_analysis_data$dge <- dge
        rnaseq_analysis_data$lrt <- lrt

        results_df <- as.data.frame(edgeR::topTags(lrt, n = Inf, sort.by = "PValue")) %>%
          tibble::rownames_to_column("Gene") %>%
          dplyr::rename(
            logFC = logFC,
            logCPM = logCPM,
            PValue = PValue,
            FDR = FDR
          ) %>%
          dplyr::select(Gene, logFC, logCPM, LR, PValue, FDR) %>%
          dplyr::arrange(FDR)

        rnaseq_analysis_data$results_table <- results_df
        rnaseq_analysis_data$status <- paste("Analyse edgeR terminée avec succès!")
        showNotification(paste("Analyse edgeR terminée!"), type = "message")

      },
      error = function(e) {
        traceback()

        rnaseq_analysis_data$dge <- NULL
        rnaseq_analysis_data$lrt <- NULL
        rnaseq_analysis_data$results_table <- NULL

        full_error_message <- as.character(if (!is.null(e$message)) e$message else "Une erreur inconnue est survenue.")
        
        display_message <- if (nchar(full_error_message) > 200) {
          paste0(substr(full_error_message, 0, 197), "...")
        } else {
          full_error_message
        }

        rnaseq_analysis_data$status <- paste("Erreur lors de l'exécution de l'analyse RNA-Seq:", display_message)

        showNotification("L'analyse RNA-Seq a échoué. Veuillez vérifier vos fichiers et la console R pour plus de détails.", type = "error", duration = 10)
      }
    )
  })

  output$rnaseq_results_table <- DT::renderDataTable({
    req(rnaseq_analysis_data$results_table)
    DT::datatable(rnaseq_analysis_data$results_table, options = list(pageLength = 10, scrollX = TRUE)) %>%
      DT::formatRound(columns = c("logFC", "logCPM", "PValue", "FDR"), digits = 3)
  })

  output$download_rnaseq_results <- downloadHandler(
    filename = function() {
      paste("rnaseq_results_edgeR_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(rnaseq_analysis_data$results_table, file, row.names = FALSE)
    }
  )

  output$mds_plot <- plotly::renderPlotly({
    req(rnaseq_analysis_data$dge)

    samples_in_analysis <- colnames(rnaseq_analysis_data$dge$counts)

    if (is.null(samples_in_analysis) || length(samples_in_analysis) < 2) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Pas assez d'échantillons pour le MDS Plot (moins de 2).", color = "red") + ggplot2::theme_void()))
    }

    metadata_for_plot <- rnaseq_analysis_data$metadata[samples_in_analysis, , drop = FALSE]

    condition_col <- input$rnaseq_condition_col
    if (!condition_col %in% names(metadata_for_plot) || length(unique(metadata_for_plot[[condition_col]])) < 2) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste0("La colonne de condition '", condition_col, "' n'a pas au moins deux niveaux uniques après filtrage pour le MDS Plot."), color = "red") + ggplot2::theme_void()))
    }

    mds_coords <- limma::plotMDS(rnaseq_analysis_data$dge, plot = FALSE)
    mds_df <- data.frame(
      Dim1 = mds_coords$x,
      Dim2 = mds_coords$y,
      Sample = rownames(metadata_for_plot),
      Condition = metadata_for_plot[[condition_col]]
    )
    plot_title_prefix <- "MDS Plot"
    x_label <- "Dimension 1"
    y_label <- "Dimension 2"

    if (is.null(mds_df) || nrow(mds_df) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Erreur inattendue: Les données pour le MDS Plot sont vides.", color = "red") + ggplot2::theme_void()))
    }

    p <- ggplot2::ggplot(mds_df, ggplot2::aes(x = Dim1, y = Dim2, color = Condition, text = paste("Sample:", Sample, "<br>Condition:", Condition))) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste(plot_title_prefix, "(Similarité des Échantillons - edgeR )"), x = x_label, y = y_label)

    plotly::ggplotly(p, tooltip = "text")
  })


  output$ma_plot <- plotly::renderPlotly({
    req(rnaseq_analysis_data$results_table, input$logfc_threshold, input$fdr_threshold)

    ma_data <- rnaseq_analysis_data$results_table %>%
      dplyr::mutate(
        is_DE = FDR < input$fdr_threshold & abs(logFC) >= input$logfc_threshold
      )

    p <- ggplot2::ggplot(ma_data, ggplot2::aes(x = logCPM, y = logFC, color = is_DE,
                                text = paste("Gene:", Gene,
                                             "<br>logFC:", round(logFC, 2),
                                             "<br>logCPM:", round(logCPM, 2),
                                             "<br>PValue:", formatC(PValue, format = "e", digits = 2),
                                             "<br>FDR:", formatC(FDR, format = "e", digits = 2)))) +
      ggplot2::geom_point(alpha = 0.6, size = 1.5) +
      ggplot2::scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
      ggplot2::geom_hline(yintercept = c(-input$logfc_threshold, input$logfc_threshold), linetype = "dashed", color = "blue") +
      ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "MA Plot (edgeR)",
                    x = "logCPM (moyenne de comptage normalisée)",
                    y = "logFC (Log2 Fold Change)") +
      ggplot2::theme(legend.position = "none")

    plotly::ggplotly(p, tooltip = "text")
  })

  output$volcano_plot <- plotly::renderPlotly({
    req(rnaseq_analysis_data$results_table, input$logfc_threshold, input$fdr_threshold)

    volcano_data <- rnaseq_analysis_data$results_table %>%
      dplyr::mutate(
        log10P = -log10(PValue),
        is_DE = FDR < input$fdr_threshold & abs(logFC) >= input$logfc_threshold
      )

    p <- ggplot2::ggplot(volcano_data, ggplot2::aes(x = logFC, y = log10P, color = is_DE,
                                text = paste("Gene:", Gene,
                                             "<br>logFC:", round(logFC, 2),
                                             "<br>PValue:", formatC(PValue, format = "e", digits = 2),
                                             "<br>FDR:", formatC(FDR, format = "e", digits = 2)))) +
      ggplot2::geom_point(alpha = 0.6, size = 1.5) +
      ggplot2::scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
      ggplot2::geom_vline(xintercept = c(-input$logfc_threshold, input$logfc_threshold), linetype = "dashed", color = "blue") +
      ggplot2::geom_hline(yintercept = -log10(input$fdr_threshold), linetype = "dashed", color = "blue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Volcano Plot (edgeR)",
                    x = "Log2 Fold Change", y = "-log10(P-value)") +
      ggplot2::theme(legend.position = "none")

    plotly::ggplotly(p, tooltip = "text")
  })

  output$heatmap_plot <- plotly::renderPlotly({
    req(rnaseq_analysis_data$dge, input$rnaseq_condition_col)

    data_matrix <- NULL
    metadata_for_heatmap <- NULL
    condition_col <- input$rnaseq_condition_col

    samples_in_analysis <- colnames(rnaseq_analysis_data$dge$counts)
    metadata_full <- rnaseq_analysis_data$metadata
    valid_samples <- intersect(samples_in_analysis, rownames(metadata_full))
    if (length(valid_samples) < 2) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Pas assez d'échantillons valides pour la Heatmap (edgeR).", color = "red") + ggplot2::theme_void()))
    }
    data_matrix <- edgeR::cpm(rnaseq_analysis_data$dge[, valid_samples], normalized.lib.sizes = TRUE, log = TRUE, prior.count = 3)
    metadata_for_heatmap <- metadata_full[valid_samples, , drop = FALSE]

    if (is.null(data_matrix) || ncol(data_matrix) < 2) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données insuffisantes pour la Heatmap (moins de 2 échantillons).", color = "red") + ggplot2::theme_void()))
    }

    sample_dist <- dist(t(data_matrix))
    sample_hclust <- hclust(sample_dist)

    ordered_samples <- sample_hclust$labels[sample_hclust$order]
    metadata_for_heatmap_ordered <- metadata_for_heatmap[ordered_samples, , drop = FALSE]
    
    if (!condition_col %in% names(metadata_for_heatmap_ordered)) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste0("La colonne de condition '", condition_col, "' n'a pas été trouvée dans les métadonnées pour la Heatmap."), color = "red") + ggplot2::theme_void()))
    }
    metadata_for_heatmap_ordered[[condition_col]] <- factor(metadata_for_heatmap_ordered[[condition_col]])

    dist_matrix <- as.matrix(sample_dist)
    dist_matrix_ordered <- dist_matrix[ordered_samples, ordered_samples]

    heatmap_df <- as.data.frame(dist_matrix_ordered) %>%
      tibble::rownames_to_column("Sample1") %>%
      tidyr::pivot_longer(
        cols = -Sample1,
        names_to = "Sample2",
        values_to = "Distance"
      ) %>%
      dplyr::mutate(
        Sample1 = factor(Sample1, levels = ordered_samples),
        Sample2 = factor(Sample2, levels = rev(ordered_samples))
      )

    num_conditions <- length(levels(metadata_for_heatmap_ordered[[condition_col]]))
    condition_colors <- RColorBrewer::brewer.pal(n = max(3, num_conditions), name = "Set1")
    names(condition_colors) <- levels(metadata_for_heatmap_ordered[[condition_col]])

    sample_conditions <- data.frame(
      Sample = rownames(metadata_for_heatmap_ordered),
      Condition = metadata_for_heatmap_ordered[[condition_col]]
    )

    heatmap_df <- heatmap_df %>%
      dplyr::left_join(sample_conditions %>% dplyr::rename(Sample1 = Sample, Condition1 = Condition), by = "Sample1") %>%
      dplyr::left_join(sample_conditions %>% dplyr::rename(Sample2 = Sample, Condition2 = Condition), by = "Sample2")


    p <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = Sample1, y = Sample2, fill = Distance,
                                    text = paste("Échantillon 1:", Sample1, "(", Condition1, ")",
                                                 "<br>Échantillon 2:", Sample2, "(", Condition2, ")",
                                                 "<br>Distance:", round(Distance, 2)))) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low = "white", high = "red") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste("Heatmap des Distances Échantillon-Échantillon (edgeR)"),
                    x = "", y = "", fill = "Distance") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = ggplot2::element_text(hjust = 1),
        axis.ticks = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )

    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(xaxis = list(autorange = TRUE, rangemode = "tozero"))
  })

  redcap_data_raw <- eventReactive(input$refresh_raw_data, {
    get_redcap_data()
  }, ignoreNULL = FALSE)

  output$year_filter_ui <- renderUI({
    data <- redcap_data_raw()
    req(data)
    if (is.data.frame(data) && nrow(data) == 0 && "record_id" %in% names(data) && !"collection_year" %in% names(data)) {
      return(NULL)
    }

    if ("collection_year" %in% names(data) && length(levels(data$collection_year)) > 0) {
      years_no_na <- sort(levels(data$collection_year)[levels(data$collection_year) != "Non renseigné"])
      choices <- c("Toutes les années", "Toutes les années (Excl. NA)", years_no_na, "Non renseigné")
      if ("Non renseigné" %in% choices) {
        choices <- c(choices[choices != "Non renseigné"], "Non renseigné")
      }
      selectInput("filter_year", "Filtrer par Année de Prélèvement :",
        choices = choices, selected = "Toutes les années"
      )
    } else {
      NULL
    }
  })

  output$genetic_diagnosis_filter_ui <- renderUI({
    data <- redcap_data_raw()
    req(data)
    if (is.data.frame(data) && nrow(data) == 0 && "record_id" %in% names(data) && !"genetic_diagnosis_name" %in% names(data)) {
      return(NULL)
    }

    if ("genetic_diagnosis_name" %in% names(data) && length(levels(data$genetic_diagnosis_name)) > 0) {
      diagnosis_choices <- sort(levels(data$genetic_diagnosis_name)[levels(data$genetic_diagnosis_name) != "Non renseigné"])
      selectInput("selected_diagnoses", "Sélectionner Diagnostic(s) :",
        choices = c("Tous les diagnostics", "Non renseigné", diagnosis_choices),
        selected = "Tous les diagnostics",
        multiple = TRUE,
        selectize = TRUE
      )
    } else {
      NULL
    }
  })

  output$disease_filter_ui_variants <- renderUI({
    data <- redcap_data_raw()
    req(data)
    if (is.data.frame(data) && nrow(data) == 0 && "record_id" %in% names(data) && !"genetic_diagnosis_name" %in% names(data)) {
      return(NULL)
    }

    if (!is.null(data) && "genetic_diagnosis_name" %in% names(data) && length(levels(data$genetic_diagnosis_name)) > 0) {
      diagnosis_choices <- sort(levels(data$genetic_diagnosis_name)[levels(data$genetic_diagnosis_name) != "Non renseigné"])
      selectInput("selected_disease_for_variants", "Choisir une ou plusieurs maladies :",
        choices = diagnosis_choices,
        selected = if (length(diagnosis_choices) > 0) diagnosis_choices[1] else NULL,
        multiple = TRUE,
        selectize = TRUE
      )
    } else {
      NULL
    }
  })

  output$origins_filter_ui_variants <- renderUI({
    data <- redcap_data_raw()
    req(data)
    if (is.data.frame(data) && nrow(data) == 0 && "record_id" %in% names(data) && !"origins" %in% names(data)) {
      return(NULL)
    }

    if (!is.null(data) && "origins" %in% names(data) && length(levels(data$origins)) > 0) {
      origin_choices <- sort(levels(data$origins)[levels(data$origins) != "Non renseigné"])
      selectInput("selected_origins_for_variants", "Choisir une ou plusieurs ethnies :",
        choices = origin_choices,
        selected = if (length(origin_choices) > 0) origin_choices[1] else NULL,
        multiple = TRUE,
        selectize = TRUE
      )
    } else {
      NULL
    }
  })


  filtered_redcap_data <- reactive({
    data <- redcap_data_raw()

    if (is.data.frame(data) && nrow(data) == 0 && "record_id" %in% names(data)) {
      return(data)
    }

    req(data)

    if (!is.null(input$filter_findings_ok) && "findings_ok" %in% names(data)) {
      if (input$filter_findings_ok == "Tous (Excl. NA)") {
        data <- data %>% dplyr::filter(findings_ok != "Non renseigné")
      } else if (input$filter_findings_ok != "Tous") {
        data <- data %>% dplyr::filter(findings_ok == input$filter_findings_ok)
      }
    }

    if (!is.null(input$filter_hospit) && "dna_during_hospit" %in% names(data)) {
      if (input$filter_hospit == "Tous (Excl. NA)") {
        data <- data %>% dplyr::filter(dna_during_hospit != "Non renseigné")
      } else if (input$filter_hospit != "Tous") {
        data <- data %>% dplyr::filter(dna_during_hospit == input$filter_hospit)
      }
    }

    if (!is.null(input$filter_year) && "collection_year" %in% names(data)) {
      if (input$filter_year == "Toutes les années (Excl. NA)") {
        data <- data %>% dplyr::filter(collection_year != "Non renseigné")
      } else if (input$filter_year != "Toutes les années") {
        data <- data %>% dplyr::filter(collection_year == input$filter_year)
      }
    }

    if (!is.null(input$selected_diagnoses) && "genetic_diagnosis_name" %in% names(data)) {
      if (!("Tous les diagnostics" %in% input$selected_diagnoses)) {
        data <- data %>% dplyr::filter(genetic_diagnosis_name %in% input$selected_diagnoses)
      }
    }

    return(data)
  })

  observe({
    data <- redcap_data_raw()
  })

  output$total_patients_box <- renderValueBox({
    data <- filtered_redcap_data()
    req(data)
    if (!"patient_unique_id" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Total Patients", icon = icon("users"), color = "red")
    } else {
      total_patients <- length(unique(data$patient_unique_id))
      valueBox(total_patients, "Total Patients", icon = icon("users"), color = "aqua")
    }
  })

  output$patients_without_ipp_box <- renderValueBox({
    data <- filtered_redcap_data()
    req(data)
    if (!"ipp_exome" %in% names(data) || !"patient_unique_id" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Patients sans IPP Renseigné", icon = icon("user-slash"), color = "red")
    } else {
        patients_without_ipp <- data %>%
            dplyr::filter(ipp_exome == "Non renseigné") %>%
            dplyr::distinct(patient_unique_id) %>%
            nrow()
        
        valueBox(patients_without_ipp, "Patients sans IPP Renseigné", icon = icon("user-slash"), color = "orange")
    }
  })

  output$avg_age_box <- renderValueBox({
    data <- filtered_redcap_data()
    req(data)
    age_data <- data %>% dplyr::filter(!is.na(age_exome), is.numeric(age_exome))
    if (!("age_exome" %in% names(data)) || all(is.na(data$age_exome)) || !is.numeric(data$age_exome) || nrow(data) == 0) {
      valueBox("N/A", "Âge Moyen (non NA)", icon = icon("calendar-alt"), color = "red")
    } else {
        unique_patients_age <- data %>%
            dplyr::filter(!is.na(age_exome), is.numeric(age_exome)) %>%
            dplyr::distinct(patient_unique_id, .keep_all = TRUE)
            
        if (nrow(unique_patients_age) == 0) {
            valueBox("N/A", "Âge Moyen (non NA)", icon = icon("calendar-alt"), color = "red")
        } else {
            avg_age <- round(mean(unique_patients_age$age_exome, na.rm = TRUE), 1)
            valueBox(avg_age, "Âge Moyen", icon = icon("calendar-alt"), color = "green")
        }
    }
  })

  output$positive_results_perc_box <- renderValueBox({
    data <- filtered_redcap_data()
    req(data)
    if (!"findings_ok" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Exomes Positifs", icon = icon("dna"), color = "red")
    } else {
      summary_data <- data %>%
        dplyr::group_by(findings_ok) %>%
        dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop')
      
      total_valid_results_patients <- summary_data %>%
        dplyr::filter(findings_ok %in% c("Yes", "No")) %>%
        dplyr::pull(count) %>%
        sum()

      positive_results_patients <- summary_data %>%
        dplyr::filter(findings_ok == "Yes") %>%
        dplyr::pull(count) %>%
        sum()

      if (total_valid_results_patients > 0) {
        percentage <- round((positive_results_patients / total_valid_results_patients) * 100, 1)
        valueBox(paste0(percentage, "%"), "Exomes Positifs", icon = icon("dna"), color = "purple")
      } else {
        valueBox("0%", "Exomes Positifs", icon = icon("dna"), color = "red")
      }
    }
  })

  output$total_prescribers_box <- renderValueBox({
    data <- filtered_redcap_data()
    req(data)
    if (!"prescriber_normalized" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Total Médecins Responsables", icon = icon("user-md"), color = "red")
    } else {
      unique_prescribers <- data %>%
        dplyr::filter(prescriber_normalized != "Non renseigné") %>%
        dplyr::distinct(patient_unique_id, .keep_all = TRUE) %>%
        dplyr::pull(prescriber_normalized) %>%
        unique() %>%
        length()

      valueBox(unique_prescribers, "Total Médecins Responsables", icon = icon("user-md"), color = "maroon")
    }
  })

  output$nanopore_seq_box <- renderValueBox({
    data <- filtered_redcap_data()
    req(data)
    if (!"nanopore_test" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Séquençage Nanopore", icon = icon("microchip"), color = "red")
    } else {
      nanopore_yes_count <- data %>%
        dplyr::filter(nanopore_test == "Yes") %>%
        dplyr::summarise(count = n_distinct(patient_unique_id)) %>%
        dplyr::pull(count)
      valueBox(nanopore_yes_count, "Séquençage Nanopore", icon = icon("microchip"), color = "blue")
    }
  })

  output$wgs_seq_box <- renderValueBox({
    data <- filtered_redcap_data()
    req(data)
    if (!"wgs_compl" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Whole Genome Sequencing", icon = icon("dna"), color = "red")
    } else {
      wgs_yes_count <- data %>%
        dplyr::filter(wgs_compl == "Yes") %>%
        dplyr::summarise(count = n_distinct(patient_unique_id)) %>%
        dplyr::pull(count)

      valueBox(wgs_yes_count, "Whole Genome Sequencing", icon = icon("dna"), color = "purple")
    }
  })


  output$sex_distribution_plot <- plotly::renderPlotly({
    data <- filtered_redcap_data()
    req(data)
    if (!"sex" %in% names(data) || nrow(data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données 'Sexe' non disponibles", color = "red") + ggplot2::theme_void()))
    }
    
    data_plot <- data %>%
      dplyr::group_by(sex) %>%
      dplyr::summarise(n = n_distinct(patient_unique_id), .groups = 'drop') %>%
      dplyr::filter(!is.na(sex))

    total_patients_for_sex_plot <- sum(data_plot$n)
    
    data_plot <- data_plot %>%
      dplyr::mutate(percentage = n / total_patients_for_sex_plot) %>%
      dplyr::arrange(sex)

    if (nrow(data_plot) == 0 || (nrow(data_plot) == 1 && data_plot$sex[1] == "Non renseigné" && data_plot$n[1] == 0)) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Aucune donnée valide pour la distribution par sexe.", color = "red") + ggplot2::theme_void()))
    }

    colors_sex <- c("Homme" = "#E47F7F", "Femme" = "#7FDBDB", "Non renseigné" = "gray")

    p <- ggplot2::ggplot(data_plot, ggplot2::aes(
      x = sex, y = n, fill = sex, text = paste("Sexe:", sex, "<br>Nombre:", n, "<br>Pourcentage:", scales::percent(percentage))
    )) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::scale_fill_manual(values = colors_sex) +
      ggplot2::geom_text(aes(label = scales::percent(percentage)), vjust = -0.5, size = 3) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Sexe", y = "Nombre de Patients Uniques", fill = "Sexe") +
      ggplot2::theme(legend.position = "none")

    plotly::ggplotly(p, tooltip = "text")
  })

  output$findings_ok_distribution_plot <- plotly::renderPlotly({
    data <- filtered_redcap_data()
    req(data)
    if (!"findings_ok" %in% names(data) || nrow(data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données 'Findings OK' non disponibles ou toutes NA", color = "red") + ggplot2::theme_void()))
    }

    data_plot <- data %>%
      dplyr::group_by(findings_ok) %>%
      dplyr::summarise(n = n_distinct(patient_unique_id), .groups = 'drop') %>%
      dplyr::filter(!is.na(findings_ok))

    total_patients_for_findings_plot <- sum(data_plot$n)

    data_plot <- data_plot %>%
      dplyr::mutate(percentage = n / total_patients_for_findings_plot) %>%
      dplyr::arrange(findings_ok)

    if (nrow(data_plot) == 0 || (nrow(data_plot) == 1 && data_plot$findings_ok[1] == "Non renseigné" && data_plot$n[1] == 0)) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Aucune donnée valide pour la distribution des résultats d'exomes.", color = "red") + ggplot2::theme_void()))
    }

    colors_status <- c("Yes" = "#2ECC71", "No" = "#E74C3C", "Non renseigné" = "gray")


    p <- ggplot2::ggplot(data_plot, ggplot2::aes(
      x = findings_ok, y = n, fill = findings_ok, text = paste("Statut:", findings_ok, "<br>Nombre:", n, "<br>Pourcentage:", scales::percent(percentage))
    )) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::scale_fill_manual(values = colors_status) +
      ggplot2::geom_text(aes(label = scales::percent(percentage)), vjust = -0.5, size = 3) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Résultat Exome (Findings OK)", y = "Nombre de Patients Uniques", fill = "Statut") +
      ggplot2::theme(legend.position = "none")

    plotly::ggplotly(p, tooltip = "text")
  })

  output$age_distribution_plot <- plotly::renderPlotly({
    data <- filtered_redcap_data()
    req(data)
    age_data <- data %>%
      dplyr::filter(!is.na(age_exome), is.numeric(age_exome)) %>%
      dplyr::distinct(patient_unique_id, .keep_all = TRUE)

    if (nrow(age_data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Pas de données d'âge valides pour l'histogramme", color = "red") + ggplot2::theme_void()))
    }

    p <- ggplot2::ggplot(age_data, ggplot2::aes(x = age_exome, text = paste("Âge:", age_exome))) +
      ggplot2::geom_histogram(binwidth = 5, fill = "#ADD8E6", color = "black") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Âge", y = "Fréquence", title = "Distribution de l'Âge des Patients")

    plotly::ggplotly(p, tooltip = "text")
  })

  output$origins_distribution_plot <- plotly::renderPlotly({
    data <- filtered_redcap_data()
    req(data)
    if (!"origins" %in% names(data) || nrow(data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données 'Origines' non disponibles.", color = "red") + ggplot2::theme_void()))
    }

    data_plot <- data %>%
      dplyr::group_by(origins) %>%
      dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop') %>%
      dplyr::filter(!is.na(origins))

    total_patients_for_origins_plot <- sum(data_plot$count)

    data_plot <- data_plot %>%
      dplyr::mutate(percentage = count / total_patients_for_origins_plot) %>%
      dplyr::arrange(dplyr::desc(count))

    if (nrow(data_plot) == 0 || (nrow(data_plot) == 1 && data_plot$origins[1] == "Non renseigné" && data_plot$count[1] == 0)) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Aucune donnée valide pour la distribution des origines.", color = "red") + ggplot2::theme_void()))
    }
    
    data_plot$origins <- factor(data_plot$origins, levels = unique(data_plot$origins))

    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      colors_origins_base <- grDevices::rainbow(nrow(data_plot))
    } else {
      num_unique_origins <- length(unique(data_plot$origins))
      if (num_unique_origins <= 8) {
        colors_origins_base <- RColorBrewer::brewer.pal(max(3, num_unique_origins), "Set2")
      } else {
        colors_origins_base <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(num_unique_origins)
      }
    }

    colors_mapping <- setNames(colors_origins_base[1:num_unique_origins], unique(data_plot$origins))
    if ("Non renseigné" %in% names(colors_mapping)) {
      colors_mapping["Non renseigné"] <- "gray"
    }


    p <- ggplot2::ggplot(data_plot, ggplot2::aes(
      x = reorder(origins, count), y = count, fill = origins,
      text = paste("Origine:", origins, "<br>Nombre:", count, "<br>Pourcentage:", scales::percent(percentage))
    )) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::scale_fill_manual(values = colors_mapping) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Origine", y = "Nombre de Patients Uniques", fill = "Origine", title = "Distribution des Origines des Patients") +
      ggplot2::theme(legend.position = "none")

    plotly::ggplotly(p, tooltip = "text")
  })

  output$patients_by_year_plot <- plotly::renderPlotly({
    data <- redcap_data_raw()
    req(data)

    if (!"collection_year" %in% names(data) || nrow(data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données d'année de prélèvement non disponibles.", color = "red") + ggplot2::theme_void()))
    }

    patients_per_year <- data %>%
      dplyr::filter(collection_year != "Non renseigné") %>%
      dplyr::group_by(collection_year) %>%
      dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop') %>%
      dplyr::arrange(collection_year)

    if (nrow(patients_per_year) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Pas de données d'année de prélèvement valides pour l'affichage.", color = "red") + ggplot2::theme_void()))
    }

    patients_per_year$collection_year_numeric <- as.numeric(as.character(patients_per_year$collection_year))

    p <- ggplot2::ggplot(patients_per_year, ggplot2::aes(
      x = collection_year_numeric, y = count, group = 1,
      text = paste("Année:", collection_year, "<br>Nombre de Patients:", count)
    )) +
      ggplot2::geom_line(color = "blue", size = 1) +
      ggplot2::geom_point(color = "darkblue", size = 3) +
      ggplot2::scale_x_continuous(breaks = patients_per_year$collection_year_numeric, labels = patients_per_year$collection_year) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Année de Prélèvement", y = "Nombre de Patients", title = "Évolution Annuelle du Nombre de Patients")

    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(
        xaxis = list(type = "category", autorange = TRUE),
        yaxis = list(autorange = TRUE, rangemode = "tozero")
      )
  })

  patients_dynamic_months <- reactive({
    data <- redcap_data_raw()
    log_data <- log_rcp_all_log_data()
    req(data, log_data)
    
    if (!("sample_date_collec" %in% names(data) && ("patient_unique_id" %in% names(data)))) {
      return(list(
        previous_month_data = tibble::tibble(IPP = character(), Nom = character(), Prénom = character(), `Date Prélèvement` = as.Date(character())),
        current_month_data = tibble::tibble(IPP = character(), Nom = character(), Prénom = character(), `Date Prélèvement` = as.Date(character())),
        previous_month_name = "Mois précédent",
        current_month_name = "Mois actuel"
      ))
    }

    today <- Sys.Date()
    current_month_num <- format(today, "%m")
    current_year_num <- format(today, "%Y")

    previous_month_date <- today %m-% months(1)
    previous_month_num <- format(previous_month_date, "%m")
    previous_year_num <- format(previous_month_date, "%Y")


    current_month_name_full <- format(today, "%B")
    previous_month_name_full <- format(previous_month_date, "%B")

    current_month_name_full <- paste0(toupper(substring(current_month_name_full, 1, 1)), substring(current_month_name_full, 2))
    previous_month_name_full <- paste0(toupper(substring(previous_month_name_full, 1, 1)), substring(previous_month_name_full, 2))


    previous_month_data <- data %>%
      dplyr::filter(
        !is.na(sample_date_collec),
        format(sample_date_collec, "%m") == previous_month_num,
        format(sample_date_collec, "%Y") == previous_year_num
      ) %>%
      dplyr::distinct(patient_unique_id, .keep_all = TRUE) %>%
      dplyr::select(IPP = ipp_exome, name, forname, sample_date_collec, first_name_referee_test, prescriber_normalized, findings_ok, genetic_diagnosis_name, patient_unique_id) %>%
      dplyr::arrange(dplyr::desc(sample_date_collec))

    current_month_data <- data %>%
      dplyr::filter(
        !is.na(sample_date_collec),
        format(sample_date_collec, "%m") == current_month_num,
        format(sample_date_collec, "%Y") == current_year_num
      ) %>%
      dplyr::distinct(patient_unique_id, .keep_all = TRUE) %>%
      dplyr::select(IPP = ipp_exome, name, forname, sample_date_collec, first_name_referee_test, prescriber_normalized, findings_ok, genetic_diagnosis_name, patient_unique_id) %>%
      dplyr::arrange(dplyr::desc(sample_date_collec))

    log_cols_to_join <- log_data %>%
      dplyr::select(patient_unique_id, RCP_Modifie, Date_Modification, Exome_Type_RCP_Value)

    previous_month_data <- previous_month_data %>%
      dplyr::left_join(log_cols_to_join, by = "patient_unique_id") %>%
      dplyr::mutate(
        RCP_Modifie = dplyr::if_else(is.na(RCP_Modifie), "Non", RCP_Modifie),
        Date_Modification = dplyr::if_else(is.na(Date_Modification), as.Date(NA), Date_Modification),
        Exome_Type_RCP_Value = dplyr::if_else(is.na(Exome_Type_RCP_Value), "Non renseigné", Exome_Type_RCP_Value)
      ) %>%
      dplyr::select(-patient_unique_id)

    current_month_data <- current_month_data %>%
      dplyr::left_join(log_cols_to_join, by = "patient_unique_id") %>%
      dplyr::mutate(
        RCP_Modifie = dplyr::if_else(is.na(RCP_Modifie), "Non", RCP_Modifie),
        Date_Modification = dplyr::if_else(is.na(Date_Modification), as.Date(NA), Date_Modification),
        Exome_Type_RCP_Value = dplyr::if_else(is.na(Exome_Type_RCP_Value), "Non renseigné", Exome_Type_RCP_Value)
      ) %>%
      dplyr::select(-patient_unique_id)


    names(previous_month_data) <- c(
      "IPP", "Nom", "Prénom", "Date Prélèvement", "Médecin Référent", 
      "Médecin Prescripteur", "Exome Positif (Findings OK)", "Nom de la Maladie Génétique",
      "RCP Modifié", "Date Dernière RCP", "Type Exome RCP Value"
    )
    names(current_month_data) <- c(
      "IPP", "Nom", "Prénom", "Date Prélèvement", "Médecin Référent", 
      "Médecin Prescripteur", "Exome Positif (Findings OK)", "Nom de la Maladie Génétique",
      "RCP Modifié", "Date Dernière RCP", "Type Exome RCP Value"
    )

    return(list(
      previous_month_data = previous_month_data,
      current_month_data = current_month_data,
      previous_month_name = previous_month_name_full,
      current_month_name = current_month_name_full
    ))
  })

  output$patients_previous_month_title <- renderText({
    paste0("Patients Ajoutés en ", patients_dynamic_months()$previous_month_name)
  })

  output$patients_current_month_title <- renderText({
    paste0("Patients Ajoutés en ", patients_dynamic_months()$current_month_name)
  })

  output$patients_previous_month_table <- DT::renderDataTable({
    data_prev_month <- patients_dynamic_months()$previous_month_data

    if (nrow(data_prev_month) == 0) {
      return(data.frame(message = paste0("Aucun patient ajouté en ", patients_dynamic_months()$previous_month_name, ".")))
    }

    DT::datatable(data_prev_month,
      options = list(pageLength = 5, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$patients_current_month_table <- DT::renderDataTable({
    data_current_month <- patients_dynamic_months()$current_month_data

    if (nrow(data_current_month) == 0) {
      return(data.frame(message = paste0("Aucun patient ajouté en ", patients_dynamic_months()$current_month_name, ".")))
    }

    DT::datatable(data_current_month,
      options = list(pageLength = 5, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$download_patients_previous_month <- downloadHandler(
      filename = function() { paste0("nouveaux_patients_mois_precedent_", Sys.Date(), ".csv") },
      content = function(file) { write.csv(patients_dynamic_months()$previous_month_data, file, row.names = FALSE) }
  )

  output$download_patients_current_month <- downloadHandler(
      filename = function() { paste0("nouveaux_patients_mois_actuel_", Sys.Date(), ".csv") },
      content = function(file) { write.csv(patients_dynamic_months()$current_month_data, file, row.names = FALSE) }
  )

  patients_for_custom_date_range <- reactive({
    data <- redcap_data_raw()
    req(data, input$download_patients_date_range)
    
    if (!"sample_date_collec" %in% names(data) || nrow(data) == 0) {
      return(tibble::tibble(IPP = character(), Nom = character(), Prénom = character(), `Date Prélèvement` = as.Date(character())))
    }

    start_date <- input$download_patients_date_range[1]
    end_date <- input$download_patients_date_range[2]

    filtered_data <- data %>%
      dplyr::filter(
        !is.na(sample_date_collec),
        sample_date_collec >= start_date,
        sample_date_collec <= end_date
      ) %>%
      dplyr::distinct(patient_unique_id, .keep_all = TRUE) %>%
      dplyr::select(
        IPP = ipp_exome, 
        name, 
        forname, 
        sample_date_collec, 
        first_name_referee_test, 
        prescriber_normalized, 
        findings_ok, 
        genetic_diagnosis_name,
        sex,
        age_exome,
        origins,
        hgvs,
        gene_name,
        nanopore_test,
        wgs_compl,
        rcp_completion_date,
        rcp_conclusion,
        exome_type_rcp
        ) %>%
      dplyr::arrange(dplyr::desc(sample_date_collec))
    
    names(filtered_data) <- c(
      "IPP", "Nom", "Prénom", "Date Prélèvement", "Médecin Référent", 
      "Médecin Prescripteur", "Exome Positif (Findings OK)", "Nom de la Maladie Génétique",
      "Sexe", "Age", "Origine", "Variant HGVS", "Nom du Gène",
      "Nanopore Test", "WGS Complété", "Date RCP", "Conclusion RCP", "Type Exome RCP"
      )

    return(filtered_data)
  })

  output$download_selected_date_range_patients <- downloadHandler(
    filename = function() {
      start_date <- format(input$download_patients_date_range[1], "%Y-%m-%d")
      end_date <- format(input$download_patients_date_range[2], "%Y-%m-%d")
      paste0("patients_custom_range_", start_date, "_to_", end_date, ".csv")
    },
    content = function(file) {
      write.csv(patients_for_custom_date_range(), file, row.names = FALSE)
    }
  )


  output$new_patients_previous_month_box <- renderValueBox({
    data_prev_month <- patients_dynamic_months()$previous_month_data
    count_prev_month <- nrow(data_prev_month)
    valueBox(count_prev_month, paste0("Nouveaux Patients en ", patients_dynamic_months()$previous_month_name), icon = icon("calendar-check"), color = "green")
  })

  output$new_patients_current_month_box <- renderValueBox({
    data_current_month <- patients_dynamic_months()$current_month_data
    count_current_month <- nrow(data_current_month)
    valueBox(count_current_month, paste0("Nouveaux Patients en ", patients_dynamic_months()$current_month_name), icon = icon("calendar-check"), color = "blue")
  })
  
  output$hospit_yes_box <- renderValueBox({
    data <- filtered_redcap_data()
    req(data)
    if (!"dna_during_hospit" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Prélevés en Hospi (Oui)", icon = icon("hospital"), color = "red")
    } else {
      nanopore_yes_count <- data %>%
        dplyr::filter(dna_during_hospit == "Yes") %>%
        dplyr::summarise(count = n_distinct(patient_unique_id)) %>%
        dplyr::pull(count)
      valueBox(nanopore_yes_count, "Prélevés en Hospi (Oui)", icon = icon("hospital"), color = "light-blue")
    }
  })

  output$hospit_no_box <- renderValueBox({
    data <- filtered_redcap_data()
    req(data)
    if (!"dna_during_hospit" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Prélevés Hors Hospi (Non)", icon = icon("clinic-medical"), color = "red")
    } else {
      count_no <- data %>%
        dplyr::filter(dna_during_hospit == "No") %>%
        dplyr::summarise(count = n_distinct(patient_unique_id)) %>%
        dplyr::pull(count)
      valueBox(count_no, "Prélevés Hors Hospi (Non)", icon = icon("clinic-medical"), color = "orange")
    }
  })

  output$dna_hospit_distribution_plot <- plotly::renderPlotly({
    data <- filtered_redcap_data()
    req(data)
    if (!"dna_during_hospit" %in% names(data) || nrow(data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données 'Prélèvement en Hospi' non disponibles ou toutes NA", color = "red") + ggplot2::theme_void()))
    }

    data_plot <- data %>%
      dplyr::group_by(dna_during_hospit) %>%
      dplyr::summarise(n = n_distinct(patient_unique_id), .groups = 'drop') %>%
      dplyr::filter(!is.na(dna_during_hospit))

    total_patients_for_hospit_plot <- sum(data_plot$n)

    data_plot <- data_plot %>%
      dplyr::mutate(percentage = n / total_patients_for_hospit_plot) %>%
      dplyr::arrange(dna_during_hospit)

    if (nrow(data_plot) == 0 || (nrow(data_plot) == 1 && data_plot$dna_during_hospit[1] == "Non renseigné" && data_plot$n[1] == 0)) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Aucune donnée valide pour la distribution des prélèvements en hospi.", color = "red") + ggplot2::theme_void()))
    }

    colors_hospit <- c("No" = "#E47F7F", "Yes" = "#7FDBDB", "Non renseigné" = "gray")

    p <- ggplot2::ggplot(data_plot, ggplot2::aes(
      x = dna_during_hospit, y = n, fill = dna_during_hospit,
      text = paste("Prélèvement:", dna_during_hospit, "<br>Nombre:", n, "<br>Pourcentage:", scales::percent(percentage))
    )) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::scale_fill_manual(values = colors_hospit) +
      ggplot2::geom_text(aes(label = scales::percent(percentage)), vjust = -0.5, size = 3) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Prélèvement en Hospitalisation", y = "Nombre de Patients Uniques", fill = "Type de Prélèvement") +
      ggplot2::theme(legend.position = "none")

    plotly::ggplotly(p, tooltip = "text")
  })

  output$prescriptions_by_prescriber_plot <- plotly::renderPlotly({
    data <- filtered_redcap_data()
    req(data)
    if (!"prescriber_normalized" %in% names(data) || nrow(data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données prescripteurs non disponibles", color = "red") + ggplot2::theme_void()))
    }

    prescriber_counts <- data %>%
      dplyr::filter(prescriber_normalized != "Non renseigné") %>%
      dplyr::group_by(prescriber_normalized) %>%
      dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop') %>%
      dplyr::arrange(dplyr::desc(count))

    if (nrow(prescriber_counts) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Pas de données de prescripteurs valides pour l'affichage.", color = "red") + ggplot2::theme_void()))
    }

    top_n <- 10
    prescriber_counts_top <- head(prescriber_counts, top_n)

    p <- ggplot2::ggplot(prescriber_counts_top, ggplot2::aes(
      x = reorder(prescriber_normalized, count), y = count,
      text = paste("Médecin:", prescriber_normalized, "<br>Prescriptions:", count)
    )) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Médecin Prescripteur", y = "Nombre de Patients Uniques", title = paste0("Top ", top_n, " Prescriptions par Médecin"))

    plotly::ggplotly(p, tooltip = "text")
  })

  output$findings_by_prescriber_table <- DT::renderDataTable({
    data <- filtered_redcap_data()
    req(data)
    if (!"prescriber_normalized" %in% names(data) || !"findings_ok" %in% names(data) || nrow(data) == 0) {
      return(data.frame(message = "Données prescripteurs ou résultats exome non disponibles."))
    }

    prescriber_findings <- data %>%
      dplyr::filter(prescriber_normalized != "Non renseigné") %>%
      dplyr::group_by(prescriber_normalized, findings_ok) %>%
      dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop_last') %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = findings_ok, values_from = count, values_fill = list(count = 0)) %>%
      dplyr::mutate(
        Total_Prescriptions = rowSums(dplyr::select(., -prescriber_normalized), na.rm = TRUE)
      ) %>%
      dplyr::rename(
        Positifs = Yes,
        Negatifs = No,
        Non_Renseignes_Findings = `Non renseigné`
      ) %>%
      dplyr::mutate(
        Total_Valid_Results = Positifs + Negatifs,
        Taux_Succes = dplyr::if_else(Total_Valid_Results > 0, round((Positifs / Total_Valid_Results) * 100, 1), 0)
      ) %>%
      dplyr::arrange(dplyr::desc(Total_Prescriptions))

    DT::datatable(prescriber_findings,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$top_prescribers_findings_plot <- plotly::renderPlotly({
    data <- filtered_redcap_data()
    req(data)
    if (!"prescriber_normalized" %in% names(data) || !"findings_ok" %in% names(data) || nrow(data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données prescripteurs ou résultats exome non disponibles", color = "red") + ggplot2::theme_void()))
    }

    top_n_prescribers <- data %>%
      dplyr::filter(prescriber_normalized != "Non renseigné") %>%
      dplyr::group_by(prescriber_normalized) %>%
      dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop') %>%
      dplyr::arrange(dplyr::desc(count)) %>%
      head(10) %>%
      dplyr::pull(prescriber_normalized)

    data_plot <- data %>%
      dplyr::filter(prescriber_normalized %in% top_n_prescribers) %>%
      dplyr::group_by(prescriber_normalized, findings_ok) %>%
      dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop_last') %>%
      dplyr::mutate(percentage = count / sum(count))

    data_plot$prescriber_normalized <- factor(data_plot$prescriber_normalized,
      levels = rev(top_n_prescribers)
    )

    colors_status <- c("Yes" = "#2ECC71", "No" = "#E74C3C", "Non renseigné" = "gray")


    p <- ggplot2::ggplot(data_plot, ggplot2::aes(
      x = prescriber_normalized, y = count, fill = findings_ok,
      text = paste("Médecin:", prescriber_normalized, "<br>Statut:", findings_ok, "<br>Nombre:", count, "<br>Pourcentage:", scales::percent(percentage))
    )) +
      ggplot2::geom_bar(stat = "identity", position = "stack") +
      ggplot2::scale_fill_manual(values = colors_status) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Médecin Prescripteur", y = "Nombre de Patients Uniques", fill = "Résultat Exome",
        title = "Répartition des Résultats Exome par Prescripteur (Top 10)"
      )

    plotly::ggplotly(p, tooltip = "text")
  })


  top_monthly_prescribers_referrers <- reactive({
    data <- redcap_data_raw()
    req(data)

    if (!("sample_date_collec" %in% names(data) && "patient_unique_id" %in% names(data) &&
          "prescriber_normalized" %in% names(data) && "first_name_referee_test" %in% names(data))) {
      return(list(
        current_month_prescribers = tibble::tibble(), previous_month_prescribers = tibble::tibble(),
        current_month_referrers = tibble::tibble(), previous_month_referrers = tibble::tibble(),
        current_month_name = "Mois actuel", previous_month_name = "Mois précédent"
      ))
    }

    today <- Sys.Date()
    current_month_start <- floor_date(today, "month")
    previous_month_start <- floor_date(today %m-% months(1), "month")

    current_month_name_full <- paste0(toupper(substring(format(current_month_start, "%B"), 1, 1)), substring(format(current_month_start, "%B"), 2))
    previous_month_name_full <- paste0(toupper(substring(format(previous_month_start, "%B"), 1, 1)), substring(format(previous_month_start, "%B"), 2))

    current_month_data <- data %>%
      dplyr::filter(
        !is.na(sample_date_collec),
        sample_date_collec >= current_month_start,
        sample_date_collec < (current_month_start %m+% months(1))
      ) %>%
      dplyr::distinct(patient_unique_id, .keep_all = TRUE)

    current_prescribers <- current_month_data %>%
      dplyr::filter(prescriber_normalized != "Non renseigné") %>%
      dplyr::count(prescriber_normalized, name = "count") %>%
      dplyr::arrange(dplyr::desc(count)) %>%
      head(5)

    current_referrers <- current_month_data %>%
      dplyr::filter(first_name_referee_test != "Non renseigné") %>%
      dplyr::count(first_name_referee_test, name = "count") %>%
      dplyr::arrange(dplyr::desc(count)) %>%
      head(5)

    previous_month_data <- data %>%
      dplyr::filter(
        !is.na(sample_date_collec),
        sample_date_collec >= previous_month_start,
        sample_date_collec < current_month_start
      ) %>%
      dplyr::distinct(patient_unique_id, .keep_all = TRUE)

    previous_prescribers <- previous_month_data %>%
      dplyr::filter(prescriber_normalized != "Non renseigné") %>%
      dplyr::count(prescriber_normalized, name = "count") %>%
      dplyr::arrange(dplyr::desc(count)) %>%
      head(5)

    previous_referrers <- previous_month_data %>%
      dplyr::filter(first_name_referee_test != "Non renseigné") %>%
      dplyr::count(first_name_referee_test, name = "count") %>%
      dplyr::arrange(dplyr::desc(count)) %>%
      head(5)

    list(
      current_month_prescribers = current_prescribers,
      previous_month_prescribers = previous_prescribers,
      current_month_referrers = current_referrers,
      previous_month_referrers = previous_referrers,
      current_month_name = current_month_name_full,
      previous_month_name = previous_month_name_full
    )
  })

  output$top_prescribers_current_month_title <- renderText({
    paste0("Top 5 Prescripteurs (", top_monthly_prescribers_referrers()$current_month_name, ")")
  })

  output$top_prescribers_previous_month_title <- renderText({
    paste0("Top 5 Prescripteurs (", top_monthly_prescribers_referrers()$previous_month_name, ")")
  })

  output$top_referrers_current_month_title <- renderText({
    paste0("Top 5 Référents (", top_monthly_prescribers_referrers()$current_month_name, ")")
  })

  output$top_referrers_previous_month_title <- renderText({
    paste0("Top 5 Référents (", top_monthly_prescribers_referrers()$previous_month_name, ")")
  })

  output$top_prescribers_current_month_plot <- plotly::renderPlotly({
    plot_data <- top_monthly_prescribers_referrers()$current_month_prescribers
    month_name <- top_monthly_prescribers_referrers()$current_month_name

    if (nrow(plot_data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste0("Aucune donnée pour le mois de ", month_name, "."), color = "gray") + ggplot2::theme_void()))
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(prescriber_normalized, count), y = count,
                                              text = paste("Médecin:", prescriber_normalized, "<br>Patients:", count))) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Médecin Prescripteur", y = "Nombre de Patients",
                    title = paste0("Top 5 Prescripteurs (", month_name, ")"))

    plotly::ggplotly(p, tooltip = "text")
  })

  output$top_prescribers_previous_month_plot <- plotly::renderPlotly({
    plot_data <- top_monthly_prescribers_referrers()$previous_month_prescribers
    month_name <- top_monthly_prescribers_referrers()$previous_month_name

    if (nrow(plot_data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste0("Aucune donnée pour le mois de ", month_name, "."), color = "gray") + ggplot2::theme_void()))
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(prescriber_normalized, count), y = count,
                                              text = paste("Médecin:", prescriber_normalized, "<br>Patients:", count))) +
      ggplot2::geom_bar(stat = "identity", fill = "darkblue") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Médecin Prescripteur", y = "Nombre de Patients",
                    title = paste0("Top 5 Prescripteurs (", month_name, ")"))

    plotly::ggplotly(p, tooltip = "text")
  })

  output$top_referrers_current_month_plot <- plotly::renderPlotly({
    plot_data <- top_monthly_prescribers_referrers()$current_month_referrers
    month_name <- top_monthly_prescribers_referrers()$current_month_name

    if (nrow(plot_data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste0("Aucune donnée pour le mois de ", month_name, "."), color = "gray") + ggplot2::theme_void()))
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(first_name_referee_test, count), y = count,
                                              text = paste("Médecin:", first_name_referee_test, "<br>Patients:", count))) +
      ggplot2::geom_bar(stat = "identity", fill = "darkseagreen") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Médecin Référent", y = "Nombre de Patients",
                    title = paste0("Top 5 Référents (", month_name, ")"))

    plotly::ggplotly(p, tooltip = "text")
  })

  output$top_referrers_previous_month_plot <- plotly::renderPlotly({
    plot_data <- top_monthly_prescribers_referrers()$previous_month_referrers
    month_name <- top_monthly_prescribers_referrers()$previous_month_name

    if (nrow(plot_data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste0("Aucune donnée pour le mois de ", month_name, "."), color = "gray") + ggplot2::theme_void()))
    }

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(first_name_referee_test, count), y = count,
                                              text = paste("Médecin:", first_name_referee_test, "<br>Patients:", count))) +
      ggplot2::geom_bar(stat = "identity", fill = "darkred") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Médecin Référent", y = "Nombre de Patients",
                    title = paste0("Top 5 Référents (", month_name, ")"))

    plotly::ggplotly(p, tooltip = "text")
  })


  output$referrals_by_referrer_plot <- plotly::renderPlotly({
    data <- filtered_redcap_data()
    req(data)
    if (!"first_name_referee_test" %in% names(data) || nrow(data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données médecins référents non disponibles", color = "red") + ggplot2::theme_void()))
    }

    referrer_counts <- data %>%
      dplyr::filter(first_name_referee_test != "Non renseigné") %>%
      dplyr::group_by(first_name_referee_test) %>%
      dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop') %>%
      dplyr::arrange(dplyr::desc(count))

    if (nrow(referrer_counts) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Pas de données de médecins référents valides pour l'affichage.", color = "red") + ggplot2::theme_void()))
    }

    top_n <- 10
    referrer_counts_top <- head(referrer_counts, top_n)

    p <- ggplot2::ggplot(referrer_counts_top, ggplot2::aes(
      x = reorder(first_name_referee_test, count), y = count,
      text = paste("Médecin Référent:", first_name_referee_test, "<br>Patients Adressés:", count)
    )) +
      ggplot2::geom_bar(stat = "identity", fill = "darkseagreen") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Médecin Référent", y = "Nombre de Patients Uniques", title = paste0("Top ", top_n, " Médecins Référents par Patients Adressés"))

    plotly::ggplotly(p, tooltip = "text")
  })

  output$findings_by_referrer_table <- DT::renderDataTable({
    data <- filtered_redcap_data()
    req(data)
    if (!"first_name_referee_test" %in% names(data) || !"findings_ok" %in% names(data) || nrow(data) == 0) {
      return(data.frame(message = "Données médecins référents ou résultats exome non disponibles."))
    }

    referrer_findings <- data %>%
      dplyr::filter(first_name_referee_test != "Non renseigné") %>%
      dplyr::group_by(first_name_referee_test, findings_ok) %>%
      dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop_last') %>%
      dplyr::ungroup() %>%
      tidyr::pivot_wider(names_from = findings_ok, values_from = count, values_fill = list(count = 0)) %>%
      dplyr::mutate(
        Total_Referrals = rowSums(dplyr::select(., -first_name_referee_test), na.rm = TRUE)
      ) %>%
      dplyr::rename(
        Positifs = Yes,
        Negatifs = No,
        Non_Renseignes_Findings = `Non renseigné`
      ) %>%
      dplyr::mutate(
        Total_Valid_Results = Positifs + Negatifs,
        Taux_Succes = dplyr::if_else(Total_Valid_Results > 0, round((Positifs / Total_Valid_Results) * 100, 1), 0)
      ) %>%
      dplyr::arrange(dplyr::desc(Total_Referrals))

    DT::datatable(referrer_findings,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$top_referrers_findings_plot <- plotly::renderPlotly({
    data <- filtered_redcap_data()
    req(data)
    if (!"first_name_referee_test" %in% names(data) || !"findings_ok" %in% names(data) || nrow(data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données médecins référents ou résultats exome non disponibles", color = "red") + ggplot2::theme_void()))
    }

    top_n_referrers <- data %>%
      dplyr::filter(first_name_referee_test != "Non renseigné") %>%
      dplyr::group_by(first_name_referee_test) %>%
      dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop') %>%
      dplyr::arrange(dplyr::desc(count)) %>%
      head(10) %>%
      dplyr::pull(first_name_referee_test)

    data_plot <- data %>%
      dplyr::filter(first_name_referee_test %in% top_n_referrers) %>%
      dplyr::group_by(first_name_referee_test, findings_ok) %>%
      dplyr::summarise(count = n_distinct(patient_unique_id), .groups = 'drop_last') %>%
      dplyr::mutate(percentage = count / sum(count))

    data_plot$first_name_referee_test <- factor(data_plot$first_name_referee_test,
      levels = rev(top_n_referrers)
    )

    colors_status <- c("Yes" = "#2ECC71", "No" = "#E74C3C", "Non renseigné" = "gray")


    p <- ggplot2::ggplot(data_plot, ggplot2::aes(
      x = first_name_referee_test, y = count, fill = findings_ok,
      text = paste("Médecin Référent:", first_name_referee_test, "<br>Statut:", findings_ok, "<br>Nombre:", count, "<br>Pourcentage:", scales::percent(percentage))
    )) +
      ggplot2::geom_bar(stat = "identity", position = "stack") +
      ggplot2::scale_fill_manual(values = colors_status) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Médecin Référent", y = "Nombre de Patients Uniques", fill = "Résultat Exome",
        title = "Répartition des Résultats Exome par Médecin Référent (Top 10)"
      )

    plotly::ggplotly(p, tooltip = "text")
  })

  output$patients_by_diagnosis_table <- DT::renderDataTable({
    data <- filtered_redcap_data()
    req(data)

    if (!"genetic_diagnosis_name" %in% names(data) || nrow(data) == 0) {
      return(data.frame(message = "Données de diagnostic génétique non disponibles."))
    }
    
    display_data <- data %>%
      dplyr::distinct(patient_unique_id, .keep_all = TRUE) %>%
      dplyr::select(IPP = ipp_exome, name, forname, sex, age_exome, genetic_diagnosis_name, findings_ok, collection_year)

    display_data <- display_data %>%
      dplyr::mutate(IPP = dplyr::if_else(IPP == "Non renseigné" | is.na(IPP), "Non renseigné", IPP))

    names(display_data) <- c("IPP Patient", "Nom", "Prénom", "Sexe", "Âge", "Diagnostic Génétique", "Findings OK", "Année Prélèvement")

    DT::datatable(display_data,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  output$total_variants_analyzed_box <- renderValueBox({
    data <- redcap_data_raw()
    req(data)
    if (!"hgvs" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Total Patients avec Info Variant", icon = icon("flask"), color = "red")
    } else {
      total_variants_patients <- data %>%
        dplyr::filter(hgvs != "Non renseigné") %>%
        dplyr::summarise(count = n_distinct(patient_unique_id)) %>%
        dplyr::pull(count)

      valueBox(total_variants_patients, "Total Patients avec Info Variant", icon = icon("flask"), color = "blue")
    }
  })

  output$hgvs_na_variants_box <- renderValueBox({
    data <- redcap_data_raw()
    req(data)
    if (!"hgvs" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Variants HGVS Non Renseignés", icon = icon("exclamation-circle"), color = "red")
    } else {
      na_hgvs_count <- data %>%
        dplyr::filter(hgvs == "Non renseigné") %>%
        dplyr::summarise(count = n_distinct(patient_unique_id)) %>%
        dplyr::pull(count)
      valueBox(na_hgvs_count, "Variants HGVS Non Renseignés", icon = icon("exclamation-circle"), color = "yellow")
    }
  })

  output$top_hgvs_variant_box <- renderValueBox({
    data <- redcap_data_raw()
    req(data)
    if (!"hgvs" %in% names(data) || nrow(data) == 0) {
      valueBox("N/A", "Variant HGVS le plus récurrent", icon = icon("fire"), color = "red")
    } else {
      variants_data <- data %>%
        dplyr::filter(hgvs != "Non renseigné") %>%
        dplyr::distinct(patient_unique_id, hgvs, .keep_all = TRUE)


      if (nrow(variants_data) == 0) {
        valueBox("N/A", "Variant HGVS le plus récurrent", icon = icon("fire"), color = "orange")
      } else {
        variant_counts <- variants_data %>%
          dplyr::group_by(hgvs) %>%
          dplyr::summarise(Count = n(), .groups = 'drop') %>%
          dplyr::mutate(Percentage = (Count / sum(Count)) * 100) %>%
          dplyr::arrange(dplyr::desc(Percentage)) %>%
          head(1)

        if (nrow(variant_counts) > 0) {
          top_variant <- variant_counts$hgvs[1]
          top_percentage <- round(variant_counts$Percentage[1], 1)

          variant_type_info <- ""
          if ("variant_type" %in% names(data)) {
            class_for_top_variant <- data %>%
              dplyr::filter(hgvs == top_variant, variant_type != "Non renseigné") %>%
              dplyr::pull(variant_type) %>%
              unique()
            if (length(class_for_top_variant) > 0) {
              variant_type_info <- paste0(" (Classe: ", paste(class_for_top_variant, collapse = ", "), ")")
            }
          }

          valueBox(
            paste0(top_variant, " (", top_percentage, "%)", variant_type_info),
            "Variant HGVS le plus récurrent (Global)",
            icon = icon("fire"),
            color = "teal"
          )
        } else {
          valueBox("N/A", "Variant HGVS le plus récurrent", icon = icon("fire"), color = "orange")
        }
      }
    }
  })

  output$recurring_variants_plot <- plotly::renderPlotly({
    data <- redcap_data_raw()
    req(data)

    required_cols <- c("hgvs", "g_hgvs", "p_hgvs", "c_hgvs")
    if (!all(required_cols %in% names(data)) || nrow(data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Données de variants (hgvs, g_hgvs, p_hgvs, c_hgvs) non disponibles.", color = "red") + ggplot2::theme_void()))
    }

    variants_data <- data %>%
      dplyr::filter(hgvs != "Non renseigné") %>%
      dplyr::distinct(patient_unique_id, hgvs, .keep_all = TRUE) %>%
      dplyr::select(hgvs, g_hgvs, p_hgvs, c_hgvs, variant_type)

    if (nrow(variants_data) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Pas de variants HGVS renseignés pour l'analyse de récurrence.", color = "red") + ggplot2::theme_void()))
    }

    variant_counts <- variants_data %>%
      dplyr::group_by(hgvs) %>%
      dplyr::summarise(
        Count = dplyr::n(),
        g_hgvs_details = paste(unique(g_hgvs[g_hgvs != "Non renseigné"]), collapse = ", "),
        p_hgvs_details = paste(unique(p_hgvs[p_hgvs != "Non renseigné"]), collapse = ", "),
        c_hgvs_details = paste(unique(c_hgvs[c_hgvs != "Non renseigné"]), collapse = ", "),
        variant_type_details = paste(unique(variant_type[variant_type != "Non renseigné"]), collapse = ", "),
        .groups = 'drop'
      ) %>%
      dplyr::mutate(
        Total_Valid_Variants = sum(Count),
        Percentage = (Count / Total_Valid_Variants) * 100
      ) %>%
      dplyr::arrange(dplyr::desc(Percentage))

    plot_data <- variant_counts %>%
      dplyr::filter(Percentage >= 1)

    if (nrow(plot_data) == 0 && nrow(variant_counts) > 0) {
      plot_data <- head(variant_counts, 5)
      showNotification("Aucun variant n'atteint le seuil de 1%. Affichage des 5 variants les plus fréquents.", type = "info", duration = 5)
    } else if (nrow(plot_data) == 0 && nrow(variant_counts) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Pas de variants HGVS renseignés pour l'analyse de récurrence.", color = "red") + ggplot2::theme_void()))
    }

    plot_data <- plot_data %>%
      dplyr::mutate(tooltip_text = paste(
        "Variant HGVS:", hgvs,
        "<br>Classe:", dplyr::if_else(variant_type_details != "", variant_type_details, "Non renseigné"),
        "<br>Occurrences:", Count,
        "<br>Pourcentage:", round(Percentage, 2), "%",
        dplyr::if_else(g_hgvs_details != "", paste0("<br>Génomique:", g_hgvs_details), ""),
        dplyr::if_else(p_hgvs_details != "", paste0("<br>Protéine:", p_hgvs_details), ""),
        dplyr::if_else(c_hgvs_details != "", paste0("<br>Codant:", c_hgvs_details), "")
      ))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = reorder(hgvs, Percentage), y = Percentage, fill = hgvs, text = tooltip_text)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Variant HGVS", y = "Pourcentage de Récurrence", title = "Variants HGVS les Plus Fréquents (Global)") +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1))

    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(yaxis = list(autorange = TRUE, rangemode = "tozero"))
  })

  output$variants_by_disease_table <- DT::renderDataTable({
    data <- redcap_data_raw()
    req(data)

    required_cols <- c("genetic_diagnosis_name", "hgvs", "g_hgvs", "p_hgvs", "c_hgvs", "variant_type")
    if (!all(required_cols %in% names(data)) || nrow(data) == 0) {
      return(data.frame(message = "Données de diagnostic ou de variants non disponibles (genetic_diagnosis_name, hgvs, g_hgvs, p_hgvs, c_hgvs, variant_type)."))
    }

    if (is.null(input$selected_disease_for_variants) || length(input$selected_disease_for_variants) == 0) {
      return(data.frame(message = "Veuillez sélectionner au moins une maladie pour afficher les variants."))
    }

    filtered_data <- data %>%
      dplyr::filter(genetic_diagnosis_name %in% input$selected_disease_for_variants) %>%
      dplyr::filter(hgvs != "Non renseigné") %>%
      dplyr::distinct(patient_unique_id, hgvs, .keep_all = TRUE)

    if (nrow(filtered_data) == 0) {
      return(data.frame(message = "Aucun variant renseigné pour la/les maladie(s) sélectionnée(s)."))
    }

    variants_per_disease <- filtered_data %>%
      dplyr::group_by(genetic_diagnosis_name, hgvs, g_hgvs, p_hgvs, c_hgvs, variant_type) %>%
      dplyr::summarise(Nombre_Occurrences = dplyr::n(), .groups = 'drop') %>%
      dplyr::arrange(genetic_diagnosis_name, dplyr::desc(Nombre_Occurrences))

    names(variants_per_disease) <- c("Maladie", "Variant HGVS", "Variant Génomique (g.hgvs)", "Variant Protéique (p.hgvs)", "Variant Codant (c.hgvs)", "Classe de Variant", "Nombre d'Occurrences")

    DT::datatable(variants_per_disease,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  get_top_single_variant_by_ethnicity <- function(data, variant_col, selected_origins) {
    req(data, variant_col, selected_origins)

    required_cols <- c(variant_col, "gene_name", "origins")
    if (!all(required_cols %in% names(data)) || nrow(data) == 0) {
      return(tibble::tibble(origins = character(), top_variant = character(), gene_name = character(), count = numeric(), percentage = numeric()))
    }

    filtered_data <- data %>%
      dplyr::filter(
        origins %in% selected_origins,
        .data[[variant_col]] != "Non renseigné",
        !is.na(.data[[variant_col]])
      ) %>%
      dplyr::distinct(patient_unique_id, .data[[variant_col]], .keep_all = TRUE)

    if (nrow(filtered_data) == 0) {
      return(tibble::tibble(origins = character(), top_variant = character(), gene_name = character(), count = numeric(), percentage = numeric()))
    }

    top_variants <- filtered_data %>%
      dplyr::group_by(origins, variant = .data[[variant_col]], gene_name) %>%
      dplyr::summarise(count = dplyr::n(), .groups = 'drop_last') %>%
      dplyr::slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
      dplyr::mutate(percentage = (count / sum(count)) * 100) %>%
      dplyr::ungroup() %>%
      dplyr::rename(top_variant = variant)

    return(top_variants)
  }

  create_top_variant_plot <- function(data_for_plot, plot_title, variant_type_label) {
    if (nrow(data_for_plot) == 0) {
      return(plotly::ggplotly(ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste("Aucune donnée pour les variants", variant_type_label, "dans les ethnies sélectionnées."), size = 4) + ggplot2::theme_void()))
    }

    p <- ggplot2::ggplot(data_for_plot, ggplot2::aes(
      x = reorder(origins, percentage), y = percentage, fill = top_variant,
      text = paste(
        "Ethnie:", origins, "<br>",
        "Top Variant (", variant_type_label, "):", top_variant, "<br>",
        "Gène:", gene_name, "<br>",
        "Pourcentage dans l'ethnie:", round(percentage, 2), "%", "<br>",
        "Occurrences:", count
      )
    )) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::geom_text(aes(label = top_variant), hjust = 0, size = 3, color = "black") +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Ethnie", y = paste("Pourcentage de Récurrence (", variant_type_label, ")", sep = ""),
        title = plot_title,
        fill = paste("Top Variant", variant_type_label)
      ) +
      ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.y = ggplot2::element_text(angle = 0, hjust = 1),
        axis.title.y = ggplot2::element_blank()
      )

    plotly::ggplotly(p, tooltip = "text") %>%
      plotly::layout(xaxis = list(autorange = TRUE, rangemode = "tozero"))
  }


  output$top_p_hgvs_by_ethnie_plot <- plotly::renderPlotly({
    data <- redcap_data_raw()
    req(data, input$selected_origins_for_variants)

    top_p <- get_top_single_variant_by_ethnicity(data, "p_hgvs", input$selected_origins_for_variants)
    create_top_variant_plot(top_p, "Top Variant Protéique (p_hgvs) par Ethnie", "p_hgvs")
  })

  output$top_c_hgvs_by_ethnie_plot <- plotly::renderPlotly({
    data <- redcap_data_raw()
    req(data, input$selected_origins_for_variants)

    top_c <- get_top_single_variant_by_ethnicity(data, "c_hgvs", input$selected_origins_for_variants)
    create_top_variant_plot(top_c, "Top Variant Codant (c_hgvs) par Ethnie", "c_hgvs")
  })

  output$top_g_hgvs_by_ethnie_plot <- plotly::renderPlotly({
    data <- redcap_data_raw()
    req(data, input$selected_origins_for_variants)

    top_g <- get_top_single_variant_by_ethnicity(data, "g_hgvs", input$selected_origins_for_variants)
    create_top_variant_plot(top_g, "Top Variant Génomique (g_hgvs) par Ethnie", "g_hgvs")
  })


  output$raw_data_table <- DT::renderDataTable({
    data <- filtered_redcap_data()
    req(data)
    if (nrow(data) == 0) {
      return(data.frame(message = "Données non disponibles ou vide."))
    }
    
    id_col_to_display <- "patient_unique_id"

    cols_to_keep <- c(
      id_col_to_display,
      "ipp_exome",
      "record_id",
      "name",
      "forname",
      "dob_exome",
      "age_exome",
      "sex",
      "origins", 
      "findings_ok",
      "dna_during_hospit",
      "sample_date_collec",
      "collection_year",
      "collection_month",
      "emaillm",
      "prescriber_normalized",
      "first_name_referee_test",
      "genetic_diagnosis_name",
      "hgvs",
      "g_hgvs",
      "p_hgvs",
      "c_hgvs",
      "variant_type",
      "gene_name",
      "nanopore_test",
      "wgs_compl"
    )

    cols_to_display_existing <- cols_to_keep[cols_to_keep %in% names(data)]
    display_data <- data[, cols_to_display_existing, drop = FALSE]

    display_names <- c()
    for (col_name in cols_to_display_existing) {
      if (col_name == "patient_unique_id") display_names <- c(display_names, "ID Patient Unique")
      else if (col_name == "ipp_exome") display_names <- c(display_names, "IPP Exome")
      else if (col_name == "record_id") display_names <- c(display_names, "Record ID REDCap")
      else if (col_name == "name") display_names <- c(display_names, "Nom")
      else if (col_name == "forname") display_names <- c(display_names, "Prénom")
      else if (col_name == "dob_exome") display_names <- c(display_names, "Date de Naissance")
      else if (col_name == "age_exome") display_names <- c(display_names, "Âge Calculé")
      else if (col_name == "sex") display_names <- c(display_names, "Sexe")
      else if (col_name == "origins") display_names <- c(display_names, "Origine")
      else if (col_name == "findings_ok") display_names <- c(display_names, "Findings OK")
      else if (col_name == "dna_during_hospit") display_names <- c(display_names, "Prélèvement Hospi")
      else if (col_name == "sample_date_collec") display_names <- c(display_names, "Date Prélèvement")
      else if (col_name == "collection_year") display_names <- c(display_names, "Année Prélèvement")
      else if (col_name == "collection_month") display_names <- c(display_names, "Mois Prélèvement")
      else if (col_name == "emaillm") display_names <- c(display_names, "Email Prescripteur (Brut)")
      else if (col_name == "prescriber_normalized") display_names <- c(display_names, "Prescripteur (Nettoyé)")
      else if (col_name == "first_name_referee_test") display_names <- c(display_names, "Référent (Homogénéisé)")
      else if (col_name == "genetic_diagnosis_name") display_names <- c(display_names, "Diagnostic Génétique")
      else if (col_name == "hgvs") display_names <- c(display_names, "Variant HGVS")
      else if (col_name == "g_hgvs") display_names <- c(display_names, "Variant Génomique")
      else if (col_name == "p_hgvs") display_names <- c(display_names, "Variant Protéique")
      else if (col_name == "c_hgvs") display_names <- c(display_names, "Variant Codant")
      else if (col_name == "variant_type") display_names <- c(display_names, "Classe de Variant")
      else if (col_name == "gene_name") display_names <- c(display_names, "Nom du Gène")
      else if (col_name == "nanopore_test") display_names <- c(display_names, "Séquençage Nanopore")
      else if (col_name == "wgs_compl") display_names <- c(display_names, "Séquençage WGS")
      else display_names <- c(display_names, col_name)
    }
    names(display_data) <- display_names

    DT::datatable(display_data,
      options = list(pageLength = 10, scrollX = TRUE),
      rownames = FALSE
    )
  })

  log_rcp_redcap_api_url <- reactiveVal("https://redcap.krctnn.com/api/")
  log_rcp_redcap_api_token <- reactiveVal("F58E14ACDF5A118B234429DC94F320F9")

  log_rcp_all_log_data <- reactive({
    req(redcap_data_raw())
    all_patient_data <- redcap_data_raw()

    output$log_rcp_status_message <- renderUI({
        p("Chargement et filtrage des mises à jour RCP en cours... Veuillez patienter.", style = "color: blue;")
    })

    tryCatch({
        rcon <- redcapConnection(url = log_rcp_redcap_api_url(), token = log_rcp_redcap_api_token())

        log_data <- exportLogging(rcon)

        patient_metadata_for_join <- all_patient_data %>%
            dplyr::select(
                record_id,
                IPP = ipp_exome,
                Nom = name,
                Prénom = forname,
                patient_unique_id
            ) %>%
            dplyr::distinct(patient_unique_id, .keep_all = TRUE)

        vars_to_filter_rcp <- c("variant_rcp", "exome_type_rcp", "rcp_daten1", "rcp_conclusion",
                                 "lecture1_clinic_wgslr", "lecture1_clinic", "apol1_rcp_rf",
                                 "rcp_inform_to_prescr", "rendu_results_bio", "incidental_findings_rcp")
        regex_pattern_rcp <- paste(vars_to_filter_rcp, collapse = "|")

        filtered_logs <- log_data[
            grepl("Update record", log_data$action, ignore.case = TRUE) & 
            !is.na(log_data$record) & log_data$record != "" &
            grepl(regex_pattern_rcp, log_data$details, ignore.case = TRUE)
        , ]
        
        if (nrow(filtered_logs) > 0) {
            required_cols_log <- c("record", "timestamp", "details")
            if (all(required_cols_log %in% names(filtered_logs))) {
                final_log_data <- filtered_logs[, required_cols_log]
                names(final_log_data)[names(final_log_data) == "record"] <- "ID_Patient_Record"
                
                final_log_data$timestamp_full <- as.POSIXct(final_log_data$timestamp, format = "%Y-%m-%d %H:%M")
                final_log_data$Date_Modification <- as.Date(final_log_data$timestamp_full)
                final_log_data$Heure_Modification <- format(final_log_data$timestamp_full, "%H:%M:%S")

                final_log_data <- final_log_data %>%
                    dplyr::left_join(patient_metadata_for_join, by = c("ID_Patient_Record" = "record_id"))
                
                final_log_data <- final_log_data %>%
                    dplyr::group_by(patient_unique_id) %>%
                    dplyr::slice_max(order_by = timestamp_full, n = 1, with_ties = FALSE) %>%
                    dplyr::ungroup()
                
                final_log_data$RCP_Modifie <- "Oui"

                final_log_data$Details_RCP_Modifies <- sapply(final_log_data$details, function(detail_string) {
                    modified_parts <- c()
                    for (var in vars_to_filter_rcp) {
                        pattern_simple_assign <- paste0(var, "\\s*=\\s*('[^']*'|\\d+|true|false|[^,]+)")
                        pattern_old_new <- paste0(var, "\\s*\\([^)]*\\)\\s*=>\\s*('[^']*'|\\d+|true|false|[^,]+)")

                        match_old_new <- regmatches(detail_string, regexpr(pattern_old_new, detail_string, perl = TRUE))
                        match_simple <- regmatches(detail_string, regexpr(pattern_simple_assign, detail_string, perl = TRUE))

                        if (length(match_old_new) > 0) {
                            modified_parts <- c(modified_parts, match_old_new)
                        } else if (length(match_simple) > 0) {
                            modified_parts <- c(modified_parts, match_simple)
                        }
                    }
                    if (length(modified_parts) > 0) {
                        return(paste(modified_parts, collapse = "; "))
                    } else {
                        return("Non spécifié / Non pertinent")
                    }
                })

                final_log_data$Exome_Type_RCP_Value <- sapply(final_log_data$details, function(detail_string) {
                    exome_val <- NA_character_
                    match_assign <- regexpr("exome_type_rcp\\s*=\\s*('[^']*'|[^,]+)", detail_string, perl = TRUE)
                    if (match_assign != -1) {
                        extracted_part <- regmatches(detail_string, match_assign)
                        exome_val <- gsub("exome_type_rcp\\s*=\\s*", "", extracted_part, perl = TRUE)
                        exome_val <- gsub("^'|'$", "", exome_val)
                    } else {
                        match_old_new_val <- regexpr("exome_type_rcp\\s*\\([^)]*\\)\\s*=>\\s*('[^']*'|[^,]+)", detail_string, perl = TRUE)
                        if (match_old_new_val != -1) {
                            extracted_part <- regmatches(detail_string, match_old_new_val)
                            exome_val <- gsub(".*=>\\s*", "", extracted_part, perl = TRUE)
                            exome_val <- gsub("^'|'$", "", exome_val)
                        }
                    }

                    if (!is.na(exome_val)) {
                        if (exome_val == "0") {
                            return("Négatif")
                        } else if (exome_val == "1") {
                            return("Positif")
                        } else if (exome_val == "2" || exome_val == "3") {
                            return("VUS")
                        } else {
                            return(exome_val)
                        }
                    }
                    return("Non renseigné")
                })


                final_log_data$details <- NULL
                final_log_data$timestamp <- NULL
                final_log_data$timestamp_full <- NULL

                final_log_data <- final_log_data %>%
                    dplyr::select(IPP, Nom, Prénom, patient_unique_id, ID_Patient_Record, Date_Modification, Heure_Modification, RCP_Modifie, Details_RCP_Modifies, Exome_Type_RCP_Value)
                
                final_log_data <- final_log_data %>%
                    dplyr::mutate(
                        IPP = dplyr::if_else(is.na(IPP) | IPP == "", "Non renseigné", IPP),
                        Nom = dplyr::if_else(is.na(Nom) | Nom == "", "Non renseigné", Nom),
                        Prénom = dplyr::if_else(is.na(Prénom) | Prénom == "", "Non renseigné", Prénom),
                        patient_unique_id = dplyr::if_else(is.na(patient_unique_id) | patient_unique_id == "", "Non renseigné", patient_unique_id)
                    )

            } else {
                final_log_data <- data.frame(
                    IPP = character(0), Nom = character(0), Prénom = character(0), patient_unique_id = character(0),
                    ID_Patient_Record = character(0), Date_Modification = as.Date(character(0)), Heure_Modification = character(0), RCP_Modifie = character(0),
                    Details_RCP_Modifies = character(0), Exome_Type_RCP_Value = character(0)
                )
            }
        } else {
            final_log_data <- data.frame(
                IPP = character(0), Nom = character(0), Prénom = character(0), patient_unique_id = character(0),
                ID_Patient_Record = character(0), Date_Modification = as.Date(character(0)), Heure_Modification = character(0), RCP_Modifie = character(0),
                Details_RCP_Modifies = character(0), Exome_Type_RCP_Value = character(0)
            )
        }
        
        output$log_rcp_status_message <- renderUI({
            p("Journal des événements de mise à jour RCP chargé et filtré avec succès !", style = "color: green;")
        })
        return(final_log_data)

    }, error = function(e) {
        output$log_rcp_status_message <- renderUI({
            p(paste("Erreur lors de la connexion ou de l'exportation du log :", e$message), style = "color: red;")
        })
        data.frame(
            IPP = character(0), Nom = character(0), Prénom = character(0), patient_unique_id = character(0),
            ID_Patient_Record = character(0), Date_Modification = as.Date(character(0)), Heure_Modification = character(0), RCP_Modifie = character(0),
            Details_RCP_Modifies = character(0), Exome_Type_RCP_Value = character(0)
        )
    })
  })

  log_rcp_filtered_data <- reactive({
      data_to_filter <- log_rcp_all_log_data()
      req(data_to_filter, input$log_rcp_date_range)

      if (!is.null(input$log_rcp_date_range) && length(input$log_rcp_date_range) >= 2) {
          start_date <- input$log_rcp_date_range[1]
          end_date <- input$log_rcp_date_range[2]
          if ("Date_Modification" %in% names(data_to_filter)) {
              data_to_filter <- data_to_filter %>%
                  dplyr::filter(Date_Modification >= start_date & Date_Modification <= end_date)
          }
      }
      
      return(data_to_filter)
  })

  output$log_rcp_log_data_table <- renderDT({
      req(log_rcp_filtered_data())
      datatable(
          log_rcp_filtered_data(),
          options = list(
              pageLength = 10,
              scrollX = TRUE,
              dom = 'Bfrtip',
              buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
          ),
          rownames = FALSE
      )
  })

  output$log_rcp_download_data <- downloadHandler(
      filename = function() {
          paste("redcap_rcp_update_export_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
          write.csv(log_rcp_filtered_data(), file, row.names = FALSE)
      }
  )

}

shinyApp(ui = ui, server = server)
