#' EGM-generated genealogical dataset of Rio Negro (deduplicated).
#'
#' The dataset contains the deduplicated records and relational fields including
#' all members of the population of Rio Negro.
#'
#' @format A data frame with 3566 rows and 11 variables: \describe{
#'   \item{ind_id}{unique universal id}
#'   \item{DoB}{date of birth}
#'   \item{DoD}{date of death}
#'   \item{sex}{sex, numeric}
#'   \item{sex2}{sex, string}
#'   \item{alive}{status in 2015}
#'   \item{father}{father id}
#'   \item{mother}{mother id}
#'   \item{children}{children id, comma separated}
#'   \item{spouses}{spouses id, comma separated}
#'   \item{id_cuest}{id of original EGM Questionnaire}
#'   }
#' @source Data collected by author in Guatemala (2015-2016).
"final"
#'
#' EGM-generated genealogical dataset of Rio Negro (with duplicates).
#'
#' The dataset contains all the records recorded in EGM interviews in Rio Negro,
#' including duplicates. Records are grouped by the EGM questionnaires in which
#' they were originally recorded.
#'
#' @format A data frame with 5803 rows and 3 variables: \describe{
#'   \item{full_id}{unique individual id including id of EGM questionnaire and
#'   line number in 'Individuals Module'}
#'   \item{idall}{unique universal id in 'final' dataset}
#'   \item{h_id}{id of original EGM Questionnaire}
#'   }
#' @source Data collected by author in Guatemala (2015-2016).
"ind.q"
#'
#' Contextual information on EGM interviews.
#'
#' The dataset contains information on the EGM interviews carried out in Rio
#' Negro.
#'
#' @format A data frame with 109 rows and 3 variables: \describe{
#'   \item{h_id}{id of original EGM Questionnaire}
#'   \item{respondent1}{line number of main respondent in
#'   'Individuals Module' of EGM questionnaire}
#'   \item{no_u_members}{total number of records in EGM questionnaire}
#'   }
#' @source Data collected by author in Guatemala (2015-2016).
"paradata"
#'
#' Time-variant family size in Rio Negro.
#'
#' The dataset contains the number of fmaily members for Rio Negro inhabitants
#' over time (1891-2016).
#'
#' @format A data frame with 70512 rows and 13 variables: \describe{
#'   \item{ego}{unique universal id}
#'   \item{year}{year of observation}
#'   \item{parents}{parents alive}
#'   \item{children}{children alive}
#'   \item{spouses}{spouses alive}
#'   \item{siblings}{sibligs alive}
#'   \item{cousins}{cousins alive}
#'   \item{uncles}{uncles alive}
#'   \item{nephews}{nephews or nieces alive}
#'   \item{grandparents}{grandparents alive}
#'   \item{grandchildren}{grandchildren alive}
#'   \item{relatives_alive_ext}{all relatives in extended family alive}
#'   \item{relatives_alive_nuc}{all relatives in nuclear family alive}
#'   }
#' @source Data collected by author in Guatemala (2015-2016).
"yearly_nets"
