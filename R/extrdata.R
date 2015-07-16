#' A subset of expression microarray of a lung dataset.
#'
#' A dataset containing a matrix of 444 samples with 150 probes.
#'
#' @format A matrix with 150 rows and 444 samples:
#' \describe{
#'   \item{row}{Probe ID}
#'   \item{col}{Sample ID}
#' }
#' @source The data is obtained from 'ClassDiscovery' package:
#'  \url{http://bioinformatics.mdanderson.org/main/OOMPA:Overview}
"lung.dataset"






#' A subset of clinical information of a lung dataset.
#'
#' A dataset containing a data frame of 444 samples with 10 clinical variables.
#'
#' @format A data frame with 444 rows and 10 variables:
#' \describe{
#'   \item{PATIENT_ID}{Patient ID}
#'   \item{SITE}{Tumor site}
#'   \item{GENDER}{Patient gender}
#'   \item{AGE_AT_DIAGNOSIS}{Diagnosis age}
#'   \item{VITAL_STATUS}{Vital status}
#'   \item{MONTHS_TO_LAST_CONTACT_OR_DEATH}{Moths to last contact or death}
#'   \item{SMOKING_HISTORY}{Smoking history}
#'   \item{PATHOLOGIC_N_STAGE}{Pathologic N stage}
#'   \item{PATHOLOGIC_T_STAGE}{Pathologic T stage}
#'   \item{Histologic.grade}{Histologic grade}
#' }
#' @source The data is obtained from 'ClassDiscovery' package:
#'  \url{http://bioinformatics.mdanderson.org/main/OOMPA:Overview}
"lung.clinical"
