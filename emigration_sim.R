
list2env(tables, envir = environment())

start <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110)
end   <- c(18, 28, 38, 48, 58, 68, 78, 88, 98, 108, 118)
names <- c('46-50', '51-55', '56-60', '61-65', '66-70', '71-75', '76-80', '81-85', '86-90', '91-96', '96-00')

age_groups <- c('< 5',  '5-9', '10-14', '15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49', '50-54', '55-59', '60-64', '65-69', '70-74', '75-79', '80-84', '85+')

col_positions <- fwf_positions(start = start, end = end, col_names = names)
## read in rates for latest time period and use only this by default
rates <- read_fwf('england_rates.txt', skip = 3, col_positions) %>%
    mutate(age_group = age_groups) %>%
    gather(time_period, rate, -age_group) %>%
    filter(time_period == '96-00') %>%
    select(-time_period)

# age group emigration numbers
age_group <- age_group %>%
    clean_names() %>%
    mutate(gender = if_else(row_number() <= 8, 'Male', 'Female')) %>%
    slice(-c(1, 9)) %>%
    gather(time_period, n, -gender_age, -gender) %>%
    rename(age = gender_age, n_age = n) %>%
    mutate(year = str_sub(time_period, 5, 8), year = parse_integer(year), time_period = NULL)

fit.extrap <- glm(n_age ~ gender + age + year, data = age_group, family = 'poisson')

# extrapolate down from 2007 down to 2000 by interpolating with poisson regression
extrap_age <- age_group %>% filter(year == 2007)

age_group_extrap <- map_dfr(c(2000:2006, 2018:2020), ~mutate(extrap_age, year = .))

pred <- broom::augment(fit.extrap, newdata = age_group_extrap, type.predict = 'response') %>%
    select(age, gender, year, n_age = .fitted) %>%
    mutate(n_age = ceiling(n_age))

emigration <- bind_rows(pred, age_group) %>% arrange(year)


## 2017 mid year estimates in northern ireland by age
pop_data <- read_excel('MYE17_SYA.xlsx', sheet = 2) %>%
    filter(
        str_detect(area,'Northern Ireland'),
        gender %in% c('Males', 'Females'),
        year %in% 2000:2017
    ) %>%
    select(age, year, MYE, gender)

## aggregate pop data into age group
pop_data <- pop_data %>%
    mutate(age_group = case_when(
               age < 18 ~ 'Less than 18 years',
               age %in% 18:24 ~ '18-24',
               age %in% 25:34 ~ '25-34',
               age %in% 35:44 ~ '35-44',
               age %in% 45:54 ~ '45-54',
               age %in% 55:64 ~ '55-64',
               age >= 65 ~ '65 years and over'
           ),
           gender = case_when(
               gender == 'Females' ~ 'Female',
               gender == 'Males' ~ 'Male'
           )) %>%
    group_by(age_group, gender, year) %>%
    summarize(n_pop = sum(MYE))

pop_data_extrap <- pop_data %>% filter(year == 2017)
pop_data_extrap <- bind_rows(
    pop_data_extrap %>% mutate(year = 2018),
    pop_data_extrap %>% mutate(year = 2019),
    pop_data_extrap %>% mutate(year = 2020))

pop_data <- pop_data %>% bind_rows(pop_data_extrap)

emigration <- emigration %>%
    rename(age_group = age) %>%
    left_join(pop_data) %>%
    mutate(rate = n_age / n_pop)


cons <- constraints
cons <- cons[ c('age_at_hire', 'time_since_first_employment', 'duration_of_employment', 'vital_status', 'sex')]

vitalc <- cons$vital_status
vitalc <- with(vitalc, list(alive = alive_confirmed + alive_assumed, dead = dead_unknown + dead_known))
cons$vital_status <- vitalc

cons <- map(cons, ~unlist(unlist(.)))

seedt <- array(1, dim = map(cons, length), dimnames = map(cons, ~names(.)))
target.list <- 1:5
weights <- Ipfp(seedt, target.list, cons)

df <- as.data.frame.table(weights[[1]]) %>% as_tibble()

## truncate replicate sample algorithm
trs <- function(table) {
    counts <- table$Freq
    rounded <- floor(counts)
    decimal <- counts - rounded
    counts <- as.integer(counts)
    topup <- sample(length(counts), round(sum(decimal)), prob = decimal)
    counts[topup] <- counts[topup] + 1
    table$integer_counts <- counts
    table
}
df <- trs(df)

cohort <- df %>% select(-Freq, -duration_of_employment) %>%
    filter(vital_status == 'alive') %>%
    select(-vital_status) %>%
    ## assume for simplicity that lowest hiring age possible is 16 and oldest is 60
    ## additionally assume that longest time since employment possible is 40 years
    mutate(
        age_at_hire = fct_recode(age_at_hire, '16-19' = '<20', '30-60' = '30+'),
        time_since_first_employment = fct_recode(time_since_first_employment, '0-20' = '<20', '30-35' = '30+')
    ) %>%
    separate(age_at_hire, into = c('aah_lb', 'aah_ub'), sep = '-') %>%
    separate(time_since_first_employment, into = c('tsfe_lb', 'tsfe_ub'), sep = '-') %>%
    mutate_at(vars(aah_lb, aah_ub, tsfe_lb, tsfe_ub), parse_integer) %>%
    mutate(
        age_lb = tsfe_lb + aah_lb,
        age_ub = tsfe_ub + aah_ub
    ) %>%
    select(age_lb, age_ub, sex, integer_counts)



## run simulations
## aggregate over sex for simplicity
emigration <- emigration %>%
    ungroup() %>%
    group_by(age_group, year) %>%
    summarize(n_age = sum(n_age),
              n_pop = sum(n_pop)) %>%
    mutate(rate = n_age / n_pop)

emigration_extreme <- filter(emigration, age_group %in% c('Less than 18 years', '65 years and over')) %>%
    ungroup() %>%
    mutate(lb = case_when(
               age_group == 'Less than 18 years' ~ 0,
               age_group == '65 years and over' ~ 65
           ),
           ub = case_when(
               age_group == 'Less than 18 years' ~ 17,
               age_group == '65 years and over' ~ Inf
           )) %>%
    select(-age_group)
emigration_normal <- filter(emigration, !(age_group %in% c('Less than 18 years', '65 years and over'))) %>%
    ungroup() %>%
    separate(age_group, into = c('lb', 'ub'), sep = '-') %>%
    mutate_at(c('lb', 'ub'), as.integer)

emigration <- bind_rows(emigration_extreme, emigration_normal)

round_to <- function(x, base) base * round(x / base)

## assume constant death rates per age group over time
## also crudely average death rates over age groups to get emigration matched groups

rates_extreme <- filter(rates, age_group %in% c('< 5', '85+')) %>%
    mutate(lb = c(0, 85), ub = c(5, Inf)) %>%
    select(-age_group)
rates_normal <- filter(rates, !(age_group %in% c('< 5', '85+'))) %>%
    separate(age_group, into = c('lb', 'ub'), sep = '-') %>%
    mutate_at(c('lb', 'ub'), as.integer)

rates <- bind_rows(rates_extreme, rates_normal)
## aggregate cohort over sex
cohort <- cohort %>%
    group_by(age_lb, age_ub) %>%
    summarize(n = sum(n)) %>%
    ungroup() %>%
    ## round to nearest 5 to simplify age group matching
    mutate(age_lb = round_to(age_lb, 5), age_ub = round_to(age_ub, 5) - 1)



age <- function(data) {
    data %>% mutate(age_lb = age_lb + 1,
                    age_ub = age_ub + 1)
}

choose_rates <- function(age_lb, age_ub) {
    rates_to_apply <- filter(rates, age_lb <= ub, lb <= age_ub)
    ## annual rate per 1k
    rate_to_apply <- mean(rates_to_apply$rate)
    rate_to_apply / 1000
}
choose_emigration <- function(age_lb, age_ub, year_cohort) {
    emigration_to_apply <- filter(emigration, year == year_cohort, age_lb <= ub, lb <= age_ub)
    emigration_rate <- sum(emigration_to_apply$n_age) / sum(emigration_to_apply$n_pop)
    emigration_rate
}
die <- function(data) {
    ## determine_rates
    local_rates <- map2_dbl(data$age_lb, data$age_ub, ~choose_rates(.x, .y))
    counts <- map_int(data$n, ~rpois(1, local_rates  * .))
    list(total_count = sum(counts), cohort = mutate(data, n = n - counts))
}

emigrate <- function(data, year) {
    local_rates <- pmap_dbl(list(data$age_lb, data$age_ub, rep(year, nrow(data))), ~choose_emigration(..1, ..2, ..3))
    counts <- map_int(data$n, ~rpois(1, local_rates * .))
    list(total_count = sum(counts), cohort = mutate(data, n = n - counts))
}

starting_cohort <- cohort

die_first <- function(starting_cohort) {
    cohort <- starting_cohort
    years <- 2000:2020
    yearly_death_counts <- length(years)
    yearly_emigration_counts <- length(years)
    ## cohort starts at 2000
    for (i in seq_along(years)) {
        ## deaths
        death_sim <- die(cohort)
        yearly_death_counts[i] <- death_sim$total_count
        cohort <- death_sim$cohort
        ## emigration
        em_sim <- emigrate(cohort, years[i])
        cohort <- em_sim$cohort
        yearly_emigration_counts[i] <- em_sim$total_count
        ## age the cohort
        cohort <- age(cohort)
    }
    list(dead = sum(yearly_death_counts), emigrated = sum(yearly_emigration_counts))
}

emigrate_first <- function(starting_cohort) {
    cohort <- starting_cohort
    years <- 2000:2020
    yearly_death_counts <- length(years)
    yearly_emigration_counts <- length(years)
    ## cohort starts at 2000
    for (i in seq_along(years)) {
        ## emigration
        em_sim <- emigrate(cohort, years[i])
        cohort <- em_sim$cohort
        yearly_emigration_counts[i] <- em_sim$total_count
        ## deaths
        death_sim <- die(cohort)
        yearly_death_counts[i] <- death_sim$total_count
        cohort <- death_sim$cohort
        ## age the cohort
        cohort <- age(cohort)
    }
    list(dead = sum(yearly_death_counts), emigrated = sum(yearly_emigration_counts))
}


