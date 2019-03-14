# functions



#' Create age groups
#'
#' \code{age_group_custom} creates age groups with given parameters (called by multiple functions)
#' @export
age_group_custom <-
  function(x, by=5, max_age = 60, first_5_different = T,
           as_factor = F, levels_custom = NA, sepa = "-") {
    # browser()

    x <- abs( as.numeric(x) )

    if(first_5_different) {
      breaks <- c(-Inf, 5, seq(5 + by,max_age,by),Inf)
      labels <- c( paste0(0, sepa, 4), paste( seq(5,max_age - by, by), seq(by + 5 - 1,max_age, by), sep =  sepa), paste0(max_age, "+"))
    } else {
      breaks <- c(-Inf, seq(by,max_age,by),Inf)
      labels <- c( paste( seq(0,max_age - by, by), seq(by - 1,max_age, by), sep =  sepa), paste0(max_age, "+"))
    }

    out <- as.character(cut(x, breaks,labels,include.lowest = TRUE, right = F))

    if (as_factor) {
      if(all(is.na(levels_custom)))
        out <- factor(out, levels =  labels)
      else
        out <- factor(out, levels =  levels_custom)
    }

    return(out)
  }

#' Census of genealogical data
#'
#' \code{census} generates census of genealogical data (includes more variables than pseudo_census function)
#' @export
census <- function(df, from, to = NA, subset_df = T, born_in_RN = F, born_bv = F,
                   gen_age = T, women_only = F, men_only = F,
                   by = 5, update_relatives = T, show_messages = F,
                   subset_on_DoB = F, min_DoB, max_DoB) {
  # If subset_df, it returns a DF with only the relevant cases
  # Otherwise, it will create a column 'keep' to indicate which values should be kept.
  if (is.na(to)) to <- from
  # keep only cases where dates are known
  # set DoD of alive people as 2010
  if (from == to) df$census <- from
  else df$census <- paste(from,to,sep = "-")
  df$DoD_temp <- df$DoD
  df$DoD_temp[df$alive == 1] <- 2100

  if (born_in_RN) df <- dplyr::filter(df, origin_literal == 1)

  if (born_bv) df <- dplyr::filter(df, origin_departamento == "baja verapaz")

  if (women_only) df <- df %>% dplyr::filter(sex == 2)

  if (men_only) df <- df %>% dplyr::filter(sex == 1)

  if (subset_df) {
    df <- df[!is.na(df$DoB) & !is.na(df$DoD_temp),]
    DoD <- df$DoD_temp
    df$DoD_temp <- NULL
    DoB <- df$DoB
    # keep only those alive in given interval
    df <- df[DoB <= to & DoD >= from,]
  }
  else {
    DoD <- df$DoD_temp
    df$DoD_temp <- NULL
    DoB <- df$DoB
    df$keep <- NA
    df$keep <- ifelse(DoB <= to & DoD >= from, T, F)
  }
  if(subset_on_DoB) {
    # keep only those alive in given interval
    df <- df %>%
      dplyr::filter(DoB >= min_DoB & DoB <= max_DoB)
  }
  if (gen_age) {
    df$age_at_census <- round(mean(c(to,from))) - df$DoB
    df$age_gr_at_census <- age_group_custom(df$age_at_census, by = by, as_factor = T)
  }
  if (update_relatives) {
    # keep original copy
    df$children_original <- df$children
    df$father_original <- df$father
    df$mother_original <- df$mother
    df$parents_original <- df$parents
    df$cousins_original <- df$cousins
    df$grandparents_original <- df$grandparents
    df$granchildren_original <- df$grandchildren
    df$spouses_original <- df$spouses
    if (show_messages) {
      print("Columns substituted for all relative types.")
      print("Original values were saved in '[relative]_original' columns.")
      # Note: chidlren cols are equivalent since the who_to_keep criteria is limited by the people actually avaiable in
      # each year-df (which has already been subset before)
      print("For spouses: Alive during census and older than 14.")
      print("For everyone else, alive during census.")
    }
  }
  return(df)
}

#' Completion rates
#'
#' \code{completion} gets completion rates, i.e. genealogical saturation.
#' @export
completion <- function(x) {
  if (length(x) > 0) {
    #browser()
    matches <- match(x,ids) # get matches
    matches <- matches[!is.na(matches)]
    if (length(matches) > 0)  ids <- ids[-matches] # remove matching pairs
    ids <<- ids # save new ids
  }
  # return percentage matched
  return(1- length(ids)/len)
}

#' Generational depth.
#'
#' \code{gen_depth} returns individual-level generational depth.
#' @export
gen_depth <- function(who) {
  # gen depth upwards
  depth_up <- 0
  current <- who
  while(length(current) > 0) {
    current <- na.omit( c(final$father[current], final$mother[current]) )
    if (length(current) > 0) depth_up <- depth_up + 1
  }

  # downwards
  depth_down <- 0
  current <- who
  while(length(current) > 0) {
    current <- na.omit( final[current,'children'] )
    current <- na.omit( unlist(strsplit(current, ";")) )
    current <- as.numeric( current[current != "NA"] )
    if (length(current) > 0) depth_down <- depth_down + 1
  }
  data.frame(up = depth_up, down = depth_down, all = depth_up + depth_down, stringsAsFactors = F)
}

#' Get pedigree object from EGM data
#'
#' \code{get_pedigree} shows pedigree object that can be plotted.
#' @export
get_pedigree <- function(famid, ind.q, final, plot = T, highlight = NA) {

  # Subset and edit data

  idall_ind <- ind.q$idall[ind.q$h_id %in% famid]

  suppressWarnings(
    ind <-
      final %>%
      dplyr::filter(ind_id %in% idall_ind) %>%
      dplyr::select(id = ind_id, sex = sex, dadid = father, momid = mother, spouses, status = alive, id_cuest) %>%
      dplyr::mutate(
        famid = 1
        , status = plyr::mapvalues(status, from = c(1,0,3, NA), to = c(0,1,3,0))
      )
  )

  ind$affected <- 0
  if(all(!is.na(highlight))) ind$affected[highlight] <- 1

  colour <-
    unlist(
      lapply(strsplit(ind$id_cuest, ";"), function(id) {
        sort(as.numeric(str_extract(id, "^[0-9]+")))[1]
      })
    )

  ind$id_cuest <- NULL

  # Note: For simplicity, keeps only one spouse in questionnaire
  l <- ind$spouses
  l[is.na(l)] <- ""
  l <- strsplit(l, ";")

  spouses <-
    unlist(
      lapply(l, function(x) {
        x <- as.numeric(x)
        out <- x[x %in% ind$id][1]
        ifelse(!length(out), NA, out)
      })
    )

  ind$spouses <- spouses

  # remove external links

  drop_d <- !( ind$dadid %in% ind$id )
  ind$dadid[drop_d] <- NA

  drop_m <- !( ind$momid %in% ind$id )
  ind$momid[drop_m] <- NA

  drop_s <- !( ind$spouses %in% ind$id )
  ind$spouses[drop_s] <- NA

  # Get all spouse relations

  sp <-
    ind %>%
    dplyr::mutate(type = 4) %>%
    select(id, spouses, type, famid) %>%
    dplyr::mutate(famid = 1) %>%
    na.omit

  sp.sort <- t(apply(sp, 1, sort))
  sp <- sp[!duplicated(sp.sort),]

  relation <- as.matrix(sp)

  # get pedigreeList object

  ped <- pedigree(id = ind$id, dadid = ind$dadid
                  , momid = ind$momid
                  , sex = ind$sex
                  , famid = ind$famid
                  , status = ind$status
                  , affected = ind$affected
                  , relation = relation
  )

  if (plot) plot.pedigree(ped['1'], col = colour, cex = 1E-4)

  return(list(ped = ped, colour = colour))

}

#' Create kinship network
#'
#' \code{kin_as_net} transforms genealogical datasets to network format.
#' @export
kin_as_net <- function(ind.q) {
  # if (ind.q$h_id[1] == 106) browser()
  from <- to <- type <- weight <- NA
  # for weighting
  weight_kin <- 1
  weight_aff <- 1
  for (n in 1:nrow(ind.q)) {
    # if (ind.q$h_id[1] == 106 & n == 97) browser()
    # browser()
    # print(n)
    #if (n ==47) browser()
    who <- ind.q$ind_id[n]
    sp <- c(ind.q$spouse1[n], ind.q$spouse2[n], ind.q$spouse3[n], ind.q$spouse4[n], ind.q$spouse5[n], ind.q$spouse6[n], ind.q$spouse7[n], ind.q$spouse8[n])
    sp <- sp[!is.na(sp)]
    f <- ind.q$father[n]
    m <- ind.q$mother[n]
    # spouse
    if (length(sp) > 0) {
      for (sp_n in sp) {
        from <- c(from, who)
        to <- c(to, sp_n)
        type <- c(type, "aff")
        weight <- c(weight,weight_aff)
      }
    }
    # father
    if (!is.na(f)) {
      from <- c(from, who)
      to <- c(to, f)
      type <- c(type, "kin")
      weight <- c(weight,weight_kin)
    }
    # mother
    if (!is.na(m)) {
      from <- c(from, who)
      to <- c(to, m)
      type <- c(type, "kin")
      weight <- c(weight,weight_kin)
    }
  }
  from <- from[!is.na(from)]
  to <- to[!is.na(to)]
  type <- type[!is.na(type)]
  weight <- weight[!is.na(weight)]
  data.frame(from,to,type,weight,stringsAsFactors = F)
}


#' Conduct 'pseudo-census' of EGM data
#'
#' \code{pseudo_census} extracts and tallies cross-sectional population data by age and sex.
#' @export
pseudo_census <- function(year, df, cuts, labels) {

  df %>%
    dplyr::mutate(
      DoD_temp = ifelse(alive == 1, 2100, DoD)
      , age_census = year - DoB
      , ag = cut(age_census, c(cuts, Inf), include.lowest= T, right = F, labels = labels)
    ) %>%
    dplyr::filter(!is.na(DoB) & !is.na(DoD_temp) & !is.na(sex2)) %>%
    dplyr::filter(DoB <= year & DoD_temp >= year) %>%
    dplyr::count(ag, sex2) %>%
    tidyr::complete(ag, sex2, fill = list(n = 0)) %>%
    tidyr::spread(sex2, n) %>%
    select(-ag)

}

#' Personalised style for ggplot
#'
#' \code{theme_dag} Provides customised text size and family.
#' @export
theme_dag <- function (base_size = 11, base_family = "", family = "",
                       legend.position = "right") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(text = element_text(size=16, family = family),
          legend.position = legend.position)
}

#' Whipple index
#'
#' \code{whipple} computes whipple index given age parameters.
#' @export
whipple <- function(data, min_age = 23, max_age = 62, printme = T) {

  data <- data[!is.na(data)]
  tmp1 <- data[data >= min_age & data <= max_age]
  tmp2 <- tmp1 %% 5
  whipple <- (length(tmp2[tmp2 == 0]) / length(tmp1)) * 500

  if(printme) {

    print("******")
    print(whipple)
    print("******")
    print("< 105 very accurate")
    print("105 to 110 relatively accurate")
    print("110 to 125 OK")
    print("125 to 175 bad")
  } else {
    whipple
  }

}

