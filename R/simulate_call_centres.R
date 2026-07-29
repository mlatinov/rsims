#' Simulate Call Centre Performance Data
#'
#' @description
#' ## Scenario
#'
#' A telecommunications company wants to understand which factors influence
#' whether a customer's issue is successfully resolved during a support call.
#'
#' For every customer call, you observe:
#'
#' * the complexity of the customer's problem,
#' * the years of experience of the support agent,
#' * whether the customer is calling about the same problem again,
#' * the duration of the call (minutes), and
#' * whether the customer's issue was successfully resolved.
#'
#' More complex problems generally require longer calls and are more likely to
#' result in repeat calls. Experienced agents typically handle calls more
#' efficiently, reducing call duration and increasing the probability of
#' resolving the customer's issue.
#'
#' Call centres are not identical. Some centres consistently achieve higher
#' resolution rates than others because of differences in management,
#' procedures, or staff training.
#'
#' The challenge is to estimate the factors affecting call resolution while
#' accounting for differences between call centres.
#'
#' ## Causal Structure
#'
#' \preformatted{
#'
#'                Problem Complexity ─────────────► Call Duration ─────────────► Resolution
#'                       │                                 │
#'                       │                                 ▼
#'                       ├────────────► Repeat Call ─────► Resolution
#'                       │
#'                       └────────────────────────────────► Resolution
#'
#' Agent Experience ──────────────────────────────► Call Duration
#'        │
#'        └──────────────────────────────────────► Resolution
#'
#' }
#'
#' Call centres differ in their baseline probability of resolving customer
#' issues through centre-specific varying intercepts.
#'
#' ## Statistical Task
#'
#' Estimate the probability that a customer's problem is resolved while
#' accounting for problem complexity, call duration, repeat calls, agent
#' experience, and differences between call centres. Compare a standard
#' logistic regression with a hierarchical logistic regression allowing
#' centre-specific intercepts.
#'
#' ## Data Generating Model
#'
#' Repeat calls are generated according to
#'
#' \deqn{
#' Repeat_i \sim Bernoulli(p_i)
#' }
#'
#' where
#'
#' \deqn{
#' logit(p_i)
#' =
#' \alpha_R
#' +
#' \beta_{Complexity} Complexity_i
#' }
#'
#' Call duration is generated as
#'
#' \deqn{
#' CallDuration_i
#' \sim
#' Normal(\mu_i,\sigma)
#' }
#'
#' where
#'
#' \deqn{
#' \mu_i
#' =
#' \alpha_C
#' +
#' \beta_{Complexity} Complexity_i
#' +
#' \beta_{Experience} Experience_i
#' }
#'
#' Finally, problem resolution is generated as
#'
#' \deqn{
#' Resolution_{ij}
#' \sim
#' Bernoulli(\pi_{ij})
#' }
#'
#' where
#'
#' \deqn{
#' logit(\pi_{ij})
#' =
#' \alpha_j
#' +
#' \beta_{Duration} CallDuration_{ij}
#' +
#' \beta_{Complexity} Complexity_{ij}
#' +
#' \beta_{Experience} Experience_{ij}
#' +
#' \beta_{Repeat} Repeat_{ij}
#' }
#'
#' with call-centre-specific intercepts
#'
#' \deqn{
#' \alpha_j
#' \sim
#' Normal(\alpha_0,\sigma_\alpha)
#' }
#'
#' @param num_center Integer. Number of call centres to simulate.
#'
#' @param complexity_score_ranges Numeric vector of length two specifying the
#' minimum and maximum possible problem complexity scores.
#'
#' @param num_call_per_center Integer. Number of customer calls generated for
#' each call centre.
#'
#' @param mean_agent_experience Numeric. Average years of experience among
#' support agents.
#'
#' @param sd_agent_experience Numeric. Standard deviation of agent experience.
#'
#' @param baseline_prob_of_repeated_call Numeric. Population-average log-odds
#' of a customer making a repeat call.
#'
#' @param r_beta_complexity Numeric. Effect of problem complexity on the
#' log-odds of a repeat call. Larger values increase the probability that
#' customers call again.
#'
#' @param baseline_call_min Numeric. Average baseline call duration (minutes)
#' before accounting for problem complexity and agent experience.
#'
#' @param c_beta_complexity Numeric. Effect of problem complexity on expected
#' call duration.
#'
#' @param c_beta_experience Numeric. Effect of agent experience on expected
#' call duration. Negative values imply that experienced agents complete calls
#' more quickly.
#'
#' @param sigma_call_min Numeric. Residual standard deviation of call duration.
#'
#' @param baseline_resolution Numeric. Population-average log-odds that a
#' customer's issue is resolved during the call.
#'
#' @param sigma_resolution Numeric. Standard deviation of the call-centre
#' specific baseline log-odds of resolution.
#'
#' @param beta_call_min Numeric. Effect of call duration on the log-odds of
#' successfully resolving the customer's issue.
#'
#' @param beta_complexity Numeric. Direct effect of problem complexity on the
#' log-odds of problem resolution.
#'
#' @param beta_agent_experience Numeric. Direct effect of agent experience on
#' the log-odds of problem resolution.
#'
#' @param beta_repeated Numeric. Effect of a repeat call on the log-odds of
#' problem resolution.
#'
#' @return A data frame with one row per customer call containing:
#'
#' \describe{
#'   \item{call_center_id}{Unique identifier of the call centre.}
#'   \item{problem_complexity}{Complexity score of the customer's issue.}
#'   \item{agent_experience}{Years of experience of the support agent.}
#'   \item{repeat_call}{Indicator for whether the customer is calling about the same issue again (0 = no, 1 = yes).}
#'   \item{call_min}{Duration of the support call in minutes.}
#'   \item{problem_resolved}{Indicator for whether the customer's issue was successfully resolved (0 = no, 1 = yes).}
#' }
#'
#' @export
simulate_call_centres <- function(
  num_center = 20,
  num_call_per_center = 50,
  complexity_score_ranges = c(0, 5),
  mean_agent_experience = 4,
  sd_agent_experience = 1.5,
  baseline_prob_of_repeated_call = -2,
  r_beta_complexity = 0.55,
  baseline_call_min = 6,
  c_beta_complexity = 1.8,
  c_beta_experience = -0.7,
  sigma_call_min = 2.5,
  baseline_resolution = 1.2,
  sigma_resolution = 0.6,
  beta_call_min = -0.08,
  beta_complexity = -0.60,
  beta_agent_experience = 0.30,
  beta_repeated = -1.10
){
  # Simulate Call centers and calls for every centers 
  n <- num_center * num_call_per_center
  call_center_id <- rep(seq_len(num_center), each = num_call_per_center)

  # Simulate the Roots .. Complexity and Agents years of experiance 
  complexity            <- runif(n, min = complexity_score_ranges[1], max = complexity_score_ranges[2])
  agent_year_experience <- pmax(rnorm(n, mean = mean_agent_experience, sd = sd_agent_experience),0) 

  # Simulate if the call is repeated Complexity -> Repeated Call 
  repeated_call <- rbinom(n, size = 1, prob = plogis(baseline_prob_of_repeated_call + r_beta_complexity * complexity))

  # Simulate the Call Minutes depending on Agent Experiance and Complexity
  mu_call_min <- baseline_call_min + c_beta_complexity * complexity + c_beta_experience * agent_year_experience 
  call_min    <- pmax(rnorm(n, mean = mu_call_min, sd = sigma_call_min), 0)  

  # Simulate if the problem is resolved
  alpha_resolution <- baseline_resolution + sigma_resolution * rnorm(num_center, mean = 0, sd = 1)   

  # Linear Predictor 
  pi_linear_predictor <- plogis(
    alpha_resolution[call_center_id] + 
    beta_call_min   * call_min   +
    beta_complexity * complexity +
    beta_agent_experience  * agent_year_experience + 
    beta_repeated * repeated_call
  )
  # Sample from Binomial Distribution with size = 1
  resolved_problem <- rbinom(n, size = 1, prob = pi_linear_predictor)

  # Combine into simulated dataset 
  sim_data <- data.frame(
    problem_resolved = resolved_problem,
    call_min         = call_min,
    repeat_call      = repeated_call,
    agent_experiance   = agent_year_experience,
    problem_complexity = complexity,
    call_center_id     = call_center_id
  )
}

