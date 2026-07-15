#' Simulate Hierarchical Student Achievement Data
#'
#' @description
#' ## Scenario
#'
#' A national education agency wants to understand which factors influence
#' standardized test performance among students across multiple states.
#'
#' Students are nested within schools, schools are nested within districts,
#' and districts are nested within states.
#'
#' For every student, you observe:
#'
#' * socioeconomic status (SES),
#' * average daily study hours,
#' * classroom size,
#' * standardized test score,
#' * the school attended,
#' * the school district, and
#' * the state.
#'
#' Students from higher socioeconomic backgrounds tend to study more each day.
#' Larger classes may reduce academic performance because teachers have less
#' time available for individual students.
#'
#' Schools, districts, and states differ in educational resources, funding,
#' teaching quality, and curriculum, producing systematic differences in
#' average student achievement. In addition, the relationship between
#' socioeconomic status and academic performance varies across states.
#'
#' The challenge is to estimate the effects of study time, class size, and
#' socioeconomic status while accounting for the hierarchical structure of the
#' educational system.
#'
#' ## Causal Structure
#'
#' \preformatted{
#'
#' Socioeconomic Status ─────────────► Study Hours ─────────────► Test Score
#'           │                                                    ▲
#'           └────────────────────────────────────────────────────┘
#'
#' Class Size ───────────────────────────────────────────────────► Test Score
#'
#' State ─────► District ─────► School ───────────────────────────► Test Score
#'
#' }
#'
#' Schools are nested within districts, and districts are nested within
#' states. Each level contributes its own variation in student achievement.
#'
#' ## Statistical Task
#'
#' Estimate the effects of socioeconomic status, study time, and class size on
#' standardized test performance while accounting for the nested educational
#' hierarchy. Compare a standard linear regression with multilevel models
#' including varying intercepts for states, districts, and schools, and
#' state-specific socioeconomic effects.
#'
#' ## Data Generating Model
#'
#' Study hours are generated according to
#'
#' \deqn{
#' StudyHours_i
#' \sim
#' Normal(\mu_i,\sigma_{Study})
#' }
#'
#' where
#'
#' \deqn{
#' \mu_i
#' =
#' \alpha_{Study}
#' +
#' \beta_{SES}SES_i
#' }
#'
#' Test scores are generated as
#'
#' \deqn{
#' TestScore_{ijkl}
#' \sim
#' Normal(\mu_{ijkl},\sigma)
#' }
#'
#' where
#'
#' \deqn{
#' \mu_{ijkl}
#' =
#' \alpha
#' +
#' \alpha_l^{State}
#' +
#' \alpha_k^{District}
#' +
#' \alpha_j^{School}
#' +
#' \beta_{Study}StudyHours_{ijkl}
#' +
#' \beta_{Class}ClassSize_{ijkl}
#' +
#' \beta_{SES,l}SES_{ijkl}
#' }
#'
#' State, district, and school intercepts are independently generated as
#'
#' \deqn{
#' \alpha_l^{State}
#' \sim
#' Normal(0,\sigma_{State})
#' }
#'
#' \deqn{
#' \alpha_k^{District}
#' \sim
#' Normal(0,\sigma_{District})
#' }
#'
#' \deqn{
#' \alpha_j^{School}
#' \sim
#' Normal(0,\sigma_{School})
#' }
#'
#' while the socioeconomic effect varies across states:
#'
#' \deqn{
#' \beta_{SES,l}
#' \sim
#' Normal(\bar{\beta}_{SES},\sigma_{\beta})
#' }
#'
#' @param num_states Integer. Number of states to simulate.
#'
#' @param num_districts_per_state Integer. Number of school districts within
#' each state.
#'
#' @param num_schools_per_district Integer. Number of schools within each
#' district.
#'
#' @param num_of_student_per_schools Integer. Number of students generated for
#' each school.
#'
#' @param class_size_range Numeric vector specifying the minimum and maximum
#' classroom size.
#'
#' @param mean_study_hours Numeric. Average daily study hours.
#'
#' @param s_beta_ses Numeric. Effect of socioeconomic status on expected study
#' hours.
#'
#' @param baseline_test_score Numeric. Population-average standardized test
#' score.
#'
#' @param sd_test_scores_states Numeric. Standard deviation of state-specific
#' intercepts.
#'
#' @param sd_test_scores_districts Numeric. Standard deviation of district
#' intercepts.
#'
#' @param sd_test_scores_schools Numeric. Standard deviation of school
#' intercepts.
#'
#' @param beta_study_h Numeric. Effect of one additional hour of study on test
#' scores.
#'
#' @param beta_class_size Numeric. Effect of one additional student in the
#' classroom on test scores.
#'
#' @param beta_bar_ses_index Numeric. Population-average effect of
#' socioeconomic status on test scores.
#'
#' @param sd_beta_ses_index Numeric. Standard deviation of the state-specific
#' socioeconomic effects.
#'
#' @param test_score_sd Numeric. Residual standard deviation of student test
#' scores.
#'
#' @param sd_study_hours Numeric. Residual standard deviation of study hours.
#'
#' @return A data frame with one row per student containing:
#'
#' \describe{
#'   \item{state_id}{Unique identifier of the state.}
#'   \item{district_id}{Unique identifier of the school district.}
#'   \item{school_id}{Unique identifier of the school.}
#'   \item{ses_index}{Standardized socioeconomic status index.}
#'   \item{study_hours}{Average daily study hours.}
#'   \item{class_size}{Number of students in the classroom.}
#'   \item{test_score}{Student standardized test score.}
#' }
#'
#' @export
simulate_schools <- function(
  num_states = 4,
  num_districts_per_state = 6,
  num_schools_per_district = 5,
  num_of_student_per_schools = 80,
  class_size_range = c(18, 32),
  mean_study_hours = 2.5,
  sd_study_hours = 0.8,
  s_beta_ses = 0.45,
  baseline_test_score = 70,
  sd_test_scores_states = 2,
  sd_test_scores_districts = 3,
  sd_test_scores_schools = 5,
  beta_study_h = 4,
  beta_class_size = -0.30,
  beta_bar_ses_index = 6,
  sd_beta_ses_index = 1,
  test_score_sd = 6
){
  ## Simulate States, Districts Schools and Students
  n_students  <- num_of_student_per_schools
  n_schools   <- num_states * num_districts_per_state * num_schools_per_district
  n_districts <- num_states * num_districts_per_state
  n <- num_states * num_districts_per_state * num_schools_per_district * num_of_student_per_schools
  
  # per-STUDENT ids 
  school_id   <- rep(seq_len(n_schools),   each = n_students)
  district_id <- rep(seq_len(n_districts), each = num_schools_per_district * n_students)
  state_id    <- rep(seq_len(num_states),  each = num_districts_per_state * num_schools_per_district * n_students)

  # Simulate Standartized SES Index and Class sizes
  ses_index  <- rnorm(n, mean = 0, sd = 1)
  class_size <- runif(n, min = class_size_range[1], max = class_size_range[2]) 

  # Simulate Study Hours as dependet on the ses index
  mu_study_h <- mean_study_hours + s_beta_ses * ses_index 
  study_h    <- pmax(rnorm(n, mean = mu_study_h, sd = sd_study_hours), 0) 

  ## Simulate Test Scores
  # Test Score Paramters 
  alpha_states    <- sd_test_scores_states    * rnorm(n = num_states, mean = 0, sd = 1)
  alpha_districts <- sd_test_scores_districts * rnorm(n = n_districts, mean = 0, sd = 1)
  alpha_schools   <- sd_test_scores_schools   * rnorm(n = n_schools, mean = 0, sd = 1)
  beta_ses        <- beta_bar_ses_index + sd_beta_ses_index * rnorm(n = num_states, mean = 0, sd = 1)

  # Build a Linear Predictor 
  mu_i <- (
    baseline_test_score 
      + alpha_states[state_id] 
      + alpha_districts[district_id] 
      + alpha_schools[school_id]
      + beta_study_h        * study_h
      + beta_class_size     * class_size
      + beta_ses[state_id] * ses_index
    )
  # Sample Test Scores from Normal Distribution 
  test_scores <- pmin(rnorm(n = n, mean = mu_i, sd = test_score_sd), 100)

  # Combine into one Simulated Dataset 
  sim_data <- data.frame(
    study_hours = study_h,
    class_size  = class_size,
    ses_index   = ses_index,
    test_scores = test_scores,
    states      = state_id,
    districts   = district_id,
    schools     = school_id
  )
}

