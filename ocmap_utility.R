library(readr)
library(tibble)
library(stringr)
library(dplyr)


## read ocmap directory
library(tidyverse)

## read ocmap rate file
library(tidyverse)


read risksets and basic processing functions

library(tidyverse)
library(survival)

## we only ever need one year calipers since every case has at least 4 controls
file <- '/Volumes/Users/adk79/Desktop/combine_riskset/OUT/lou/Cancer/1_year_caliper/0101RSET.OUT'

full_data_file <- '/Volumes/Users/adk79/Desktop/combine_local/DAT/lou/combined.V20'
source('read_ocmap_2.R')
full_data <- read_ocmap_file(full_data_file, rt1_specd, rt2_spec, rt3_spec)
full_data_simp <- full_data %>% select(id, cod, cd_average_intensity, cd_cum_exposure) %>% distinct

old_data <- '/Volumes/Users/Projects/Chloroprene/CommonFiles/OCMAP Files/CD-L_OF7_DEIDENTIFIED.V02'
old_data <- read_ocmap_file(old_data, rt1_specd, rt2_spec, rt3_spec)
old_data <- old_data %>% mutate(job_end_date = lubridate::mdy(job_end_date))
old_data %>% filter(job_end_date == max(job_end_date))

spec <- tribble(
    ~name, ~start, ~stop,
    'caseno', 1, 6,
    'ccind', 8, 8,
    'recno', 10, 15,
    'sex', 17, 17,
    'race', 19, 19,
    'icda', 21, 24,
    'mob', 26, 27,
    'dob', 29, 30,
    'yob', 32, 35,
    'moh', 37, 38,
    'doh', 40, 41,
    'yoh', 43, 46,
    'mos', 48, 49,
    'dos', 51, 52,
    'yos', 54, 57,
    'most', 59, 60,
    'dost', 62, 63,
    'yost', 65, 68,
    'vs', 70, 71,
    ## 'plant', 72, 73,
    ## 'moev', 75, 76,
    ## 'doev', 78, 79,
    ## 'yoev', 81, 84,
    ## 'moen', 86, 87,
    ## 'doen', 89, 90,
    ## 'yoen', 92, 95,
    ## 'mosep', 97, 98,
    ## 'dosep', 100, 101,
    ## 'yosep', 103, 106,
    ## 'isdead', 108, 108,
    ## 'waswork', 110, 110,
    'dslw', 112, 117,
    'plantev', 119, 122,
    'eagedays', 124, 129,
    'latdays', 131, 136,
    'nosel', 138, 140,
    'noexp', 142, 144,
    'noavg', 146, 148,
    'plant', 150, 165,
    'doe', 166, 181,
    ## 'yoh', 182, 197,
    'aah', 198, 213,
    ## 'cdexposedever', 214, 229,
    ## 'vcexposedever', 230, 245,
    'duration_employment', 246, 261,
    'cum_cd', 263, 278,
    'cum_vc', 279, 294,
    'noCDnoVC', 295, 310,
    'CDnoVC', 311, 326,
    'noCDVC', 327, 342,
    'CDVC', 343, 358,
    'cum_CDVC', 359, 374,
    'cum_CDnoVC', 375, 390,
    'dur_cd', 391, 406,
    'dur_vc', 407, 422,
    'bc', 423, 438,
    'wc', 439, 454
    )

data <- read_fwf(file, fwf_positions(start = spec$start, end = spec$stop, col_names = spec$name))



datat <- data %>% mutate(
                      race = factor(race),
                      sex = factor(sex),
                      eventage = eagedays / 365.25,
                      latency=latdays/365.25,
                      duration_employment= duration_employment / 365.25,
                      cum_cd = cum_cd / 365.25,
                      cum_vc = cum_vc / 365.25,
                      noCDnoVC = noCDnoVC / 365.25,
                      CDnoVC = CDnoVC / 365.25,
                      ## noCDVC = noCDVC / 365.25,
                      cum_CDVC = cum_CDVC / 365.25,
                      cum_CDnoVC = cum_CDnoVC / 365.25,
                      dur_cd = dur_cd / 365.25,
                      dur_vc = dur_vc / 365.25,
                      cd_aie = if_else(dur_cd != 0, cum_cd / dur_cd, 0),
                      vc_aie = if_else(dur_vc != 0, cum_vc / dur_vc, 0),
                      bc = bc / 365.25,
                      wc = wc / 365.25,
                      white_collar = wc > bc)
                     


    
lou_data <- datat %>% filter(plant == '3')
pon_data <- datat %>% filter(plant == '1')
lou_cases <- lou_data %>% filter(ccind == 1)
pon_cases <- pon_data %>% filter(ccind == 1)
lou_data <- semi_join(lou_data, lou_cases, by = 'caseno')
pon_data <- semi_join(pon_data, pon_cases, by = 'caseno')



split_risksets <- function(data) {
    any_dupes <- data %>% filter(ccind == 1) %>%
        count(caseno) %>%
        count(n) %>%
        pull(nn) %>%
        length
    if(any_dupes == 1) {
        return(data)
    }
            

    data <- data %>% mutate(bd = lubridate::mdy(str_c(mob, dob, yob)))
    counts <- data %>% filter(ccind == 1) %>%  count(caseno)
    singles <- filter(counts, n == 1)
    multiples <- filter(counts, n > 1)
    singles <- data %>% semi_join(singles, by = 'caseno')
    multiples <- data %>% semi_join(multiples, by = 'caseno')
    multiple_cases <- multiples %>% filter(ccind == 1)
    multiple_controls <- multiples %>% filter(ccind == 2)

    new_cases <- vector(mode = 'list', length = nrow(multiple_cases))
    for (i in 1:nrow(multiple_cases)) {
        case <- multiple_cases[i,]
        eligible_controls <- multiple_controls %>% filter(caseno == case$caseno)
        controls <- eligible_controls %>% filter(abs(as.numeric(bd - case$bd)) <= 365)

        new_case <- bind_rows(case, controls) %>%
            mutate(caseno = str_c(caseno, '.', i))
        new_cases[[i]] <- new_case
    }

    new_cases <- reduce(new_cases, rbind)
    rbind(singles, new_cases)
}

lou_data <- split_risksets(lou_data)
pon_data <- split_risksets(pon_data)





dur_exposures <- c(
    'latency', 'duration_employment',
    'dur_vc', 'dur_cd',
    'bc', 'wc'
)

cat_exposures <-
    c('sex', 'race', 'white_collar')

exposures <- c(
    'cum_cd', 
    'noCDnoVC', 'CDnoVC',
    'cum_CDVC',  ##noCDVC,
    'CDVC',
    'cum_CDnoVC',
    'cd_aie'
)
exposures_vc <- c(
    'cum_vc',
    'vc_aie'
)


## logic for preparing the various exposures for analysis 
get_cuts <- function(ccind, exposure) {
    quants <- exposure[ccind == 1] %>%
        quantile(c(0.25, 0.5, 0.75))

    cut(exposure, breaks = unique(c(0, quants, Inf)), include.lowest = TRUE, right = FALSE)
}



total_exposures <- c(dur_exposures, cat_exposures, exposures)


processor <- function(data) {
    data <- data %>% mutate(sex = if_else(sex == 1, 'Male', 'Female'),
                    race = if_else(race == 1, 'White', 'NonWhite'),
                    white_collar = if_else(white_collar, 'TRUE', 'FALSE'))
    data <- data %>% mutate(latency = case_when(
                        latency < 20 ~ '<20',
                        latency < 30 ~ '20 - 29',
                        latency >= 30 ~ '30+'),
                    duration_employment = case_when(
                        duration_employment < 5 ~ '<5',
                        duration_employment < 20 ~ '5-20',
                        duration_employment >= 20 ~ '20+'),
                    dur_cd = case_when(
                        dur_cd < 10 ~ '<10',
                        dur_cd < 20 ~ '10-19',
                        dur_cd >= 20 ~ '20+'),
                    dur_vc = case_when(
                        dur_vc < 5 ~ '<5',
                        dur_vc < 10 ~ '5-10',
                        dur_vc >= 10 ~ '10+'),
                    bc = case_when(
                        bc < 5 ~ '<5',
                        bc < 10 ~ '5-10',
                        bc >= 10 ~ '10+'
                    ),
                    bc2 = case_when(
                        bc < 10 ~ '<10',
                        bc < 20 ~ '10-20',
                        bc >= 20 ~ '20+'
                    ),
                    wc = case_when(
                        bc < 5 ~ '<5',
                        bc < 10 ~ '5-10',
                        bc >= 10 ~ '10+'
                    )
                    )


    data <- data %>%
        mutate_at(vars(one_of(exposures)), function(exposure) get_cuts(data$ccind, exposure))
    data
    
}

## louisville
models <- total_exposures %>%
    map(~str_c('ccind == 1 ~ strata(caseno) + ', .)) %>%
    map(as.formula) %>%
    map(safely(function(.) clogit(., data = processor(lou_data))))

valid_results <- map_lgl(models, ~!is.null(pluck(., 'result')))

coefs <- models[valid_results] %>% map2(total_exposures[valid_results], ~mutate(broom::tidy(.x$result), exposure = .y)) %>% bind_rows() 


globals <- models[valid_results] %>%
    map(~pluck(., 'result')) %>%
    keep(~length(coef(.)) > 1) %>%
    map(anova) %>%
    map(broom::tidy) %>%
    bind_rows() %>%
    filter(!is.na(statistic))

coefs %>% knitr::kable()

globals %>% knitr::kable()


## pontchartrain




process_table <- function(table, age_groups, year_groups, sex, race, num_icd_revs = 5) {

    descript <- slice(table, 1)
    cod <- str_sub(descript$value, 30)

    ## this doesn't work generally at the moment, figure out where this needs to be pulled from
    icd_revs <- slice(table, 2:(2 + num_icd_revs))
    rates <- slice(table, (3 + num_icd_revs):(3 + num_icd_revs + length(age_groups))) %>%
        mutate(value = map_chr(value, ~str_sub(.,10,)))

    split_string <- function(string) {
        num_rates <- vector('double', length(year_groups))
        for (i in 1:length(year_groups)) {
            num_rates[i] <- str_sub(string, 10 * (i-1) + 1, 10 * (i - 1) + 1 + 9) %>% parse_double()
        }
        names(num_rates) <- year_groups
        num_rates
    }

    rates <- rates %>%
        mutate(value = map(value, split_string)) %>%
        pull(value) %>%
        do.call(rbind, .) %>%
        as_tibble() %>%
        mutate(age_groups = age_groups,
               sex = sex,
               race = race,
               cod = cod)

    list(icd_revs = mutate(icd_revs, cod = cod), rates = rates)
}


read_rate_file <- function(source_file, sex, race) { 
    rates_text <- read_lines(source_file)

    num_groups <- rates_text[[1]] %>% str_sub(1,3) %>% parse_integer()

    cl <- rates_text[[2]]
    num_entry <- str_sub(cl, 1, 3) %>% parse_integer()

    age_groups <- vector('character', num_entry)
    for (i in 1:num_entry) {
        age_groups[i] <- str_sub(cl, 5*i - 1, 5*i + 3)
    }

    cl <- rates_text[[3]]
    num_entry <- str_sub(cl, 1, 3) %>% parse_integer()
    year_groups <- vector('character', num_entry)
    for (i in 1:num_entry) {
        year_groups[i] <- str_sub(cl, 5*i - 1, 5*i + 3)
    }

    data <- enframe(rates_text[4:length(rates_text)])
    data <- data %>% mutate(group_num = str_sub(value, 1, 3)) %>%
        group_by(group_num)

    stopifnot(n_groups(data) == num_groups)

    tables <- data %>% group_split()

    processed <- tables %>%
        modify(~process_table(., age_groups, year_groups, sex, race))
    rates <- modify(processed, 'rates')
    icd_revs <- modify(processed, 'icd_revs')
    
    rates <- rates %>%
        reduce(bind_rows) %>%
        gather(year, rate, -age_groups, -sex, -race, -cod)
    icd_revs <- reduce(icd_revs, bind_rows) %>% select(-group_num, -name)
    out = list(rates = rates, icd_revs = icd_revs)
    out
}

icd_revs <-     read_rate_file('/Volumes/Users/adk79/Desktop/combine_local/DAT/lou/lou.NWF', sex = 'f', race = 'nw')$icd_revs
rates <- bind_rows(
    read_rate_file('/Volumes/Users/adk79/Desktop/combine_local/DAT/lou/lou.NWF', sex = 'f', race = 'nw')$rates,
    read_rate_file('/Volumes/Users/adk79/Desktop/combine_local/DAT/lou/lou.WF', sex = 'f', race = 'w')$rates,
    read_rate_file('/Volumes/Users/adk79/Desktop/combine_local/DAT/lou/lou.WM', sex = 'm', race = 'w')$rates,
    read_rate_file('/Volumes/Users/adk79/Desktop/combine_local/DAT/lou/lou.NWM', sex = 'm', race = 'nw')$rates,
    ## read_rate_file('/Volumes/Users/adk79/Desktop/combine_local/DAT/lou.WF', sex = 'f', race = 'w'),
    ## read_rate_file('/Volumes/Users/adk79/Desktop/combine_local/DAT/lou.WF', sex = 'f', race = 'w'),
    ## read_rate_file('/Volumes/Users/adk79/Desktop/combine_local/DAT/lou.WF', sex = 'f', race = 'w'),
    ## read_rate_file('/Volumes/Users/adk79/Desktop/combine_local/DAT/lou.WF', sex = 'f', race = 'w'),
)

rates <- rates %>%    mutate(
        age_groups = case_when(
            age_groups == "  < 5" ~ 0L,
            age_groups == "85+" ~ 85L,
            TRUE ~ parse_integer(str_sub(age_groups, 1,2))
        ),
        year = parse_integer(str_sub(year, 1, 2)) + 1900L
    )



read_directory <- function(project_name, directory) {
cohort_configs <- read_lines(str_c(directory, 'SYS/', project_name, '/', project_name, '.CID')) %>%
    enframe() %>%
    rename(cname = name, cvalue = value) %>%
    mutate(cvalue = str_remove_all(cvalue, '[^A-Za-z0-9]'), 
           cname = cname - 1,
           cname = as.character(cname),
           cname = if_else(str_length(cname) == 1, str_c('0', cname), cname))

task_configs <- read_lines(str_c(directory, 'SYS/', project_name, '/', project_name, '.TID')) %>%
    enframe() %>%
    rename(tname = name, tvalue = value) %>%
    mutate(tvalue = str_remove_all(tvalue, '[^A-Za-z0-9]'),
           tname = tname - 1,
           tname = as.character(tname),
           tname = if_else(str_length(tname) == 1, str_c('0', tname), tname))


configs <- crossing(cohort_configs, task_configs)

smr_files <- configs %>%
    mutate(file_name = str_c(cname, tname, 'RPT.TT4'),
           full_file_name = str_c(directory, 'OUT/', project_name, '/', file_name),
           file_text = map(full_file_name, safely(~read_lines(.))),
           successful = map(file_text, ~pluck(., 'error')),
           successful = map_lgl(successful, is.null)) %>%
    filter(successful) %>%
    mutate(file_text = map(file_text, ~pluck(., 'result'))) %>%
    select(-full_file_name, -successful, -cname, -tname)


process_tt4_file <- function(lines) {
    entries <- tibble()
    start <- 1
    repeat {

        repeat{
            
            if(start > length(lines)) {
                return(entries)
            }
            
            if(str_detect(lines[[start]], 'Project')) {
                break;
            }


            start <- start + 1
        }

        stop <- start
        repeat {
            if (str_detect(lines[[stop]], 'SIGNIFICANT AT 5')) {
                break;
            }
            stop <- stop + 1
        }

        block <- lines[start:stop]
        top_info <- tibble(plant = str_sub(block[1], 71, 72),
                           exp1 = str_sub(block[1], 96, 106),
                           exp2 = str_sub(block[2], 96, 106),
                           exp3 = str_sub(block[3], 96, 106),
                           exp4 = str_sub(block[4], 96, 106),
                           exp5 = str_sub(block[1], 120, 133),
                           exp6 = str_sub(block[2], 120, 133),
                           exp7 = str_sub(block[3], 120, 133),
                           exp8 = str_sub(block[4], 120, 133),
                           race = str_sub(block[2], 71, 82),
                           sex = str_sub(block[3], 71, 82),
                           latency = str_sub(block[4], 71, 82),
                           age = str_sub(block[5], 71, 82),
                           time = str_sub(block[6], 71, 82),
                           no_at_risk = str_sub(block[6], 95, 110),
                           pyrs = str_sub(block[6], 120, 133),
                           unknown = str_sub(block[length(block)], 61, 62))
        
        smr_table <- block[16:(length(block)-2)]

        ## smr_table %>% map(
        smr_table <- smr_table %>% keep(~str_detect(., '[A-Za-z]')) %>% .[1:(length(.) - 1)]

        smr_table <- smr_table %>% map_dfr(~tibble(
                                               cod = str_sub(., 1, 53),
                                               obs = str_sub(., 56, 64),
                                               exp = str_sub(., 68, 76),
                                               smr = str_sub(., 80, 87),
                                               sig = str_sub(., 89, 90),
                                               lower_95 = str_sub(., 95, 104),
                                               upper_95 = str_sub(., 111, 119),
                                               lower_99 = str_sub(., 125, 134),
                                               upper_99 = str_sub(., 141, 149)))

        out <- bind_cols(smr_table, top_info[rep(1, nrow(smr_table)) , ])
        entries <- bind_rows(entries, out)
        start <- stop + 1
    }
    
}



tabs <- smr_files %>% mutate(table = map(file_text, process_tt4_file)) %>%
    unnest(table) %>%
    mutate(obs = parse_integer(obs),
           exp = parse_double(exp),
           smr = parse_double(smr),          
           lower_95 = parse_double(lower_95),
           upper_95 = parse_double(upper_95),
           lower_99 = parse_double(lower_99),
           upper_99 = parse_double(upper_99)) %>%
    mutate_if(is.character, str_squish) %>%
    mutate_if(is.character, parse_guess)

}


## read ocmap file
## given a row specification (as a tibble)
## parse a row and turn into a row for final dataset
parse_row_with_spec <- function(spec, row) {
    spec %>%
        ## possible speedup here if str_sub is doing redundant work
        mutate(entry = str_sub(row, start, stop)) %>%
        select(name, entry) %>%
        deframe
}


read_ocmap_file <- function(file, rt1_spec = rt1_specd, rt2_spec, rt3_spec, try_convert = TRUE) {

    raw_data <- read_lines(file)

    row_number <- 1
    line_number <- 1

    ## pre initialize more entries than are necessary to speed things up
    entries <- matrix('', nrow = length(raw_data), ncol = nrow(rt1_spec) + nrow(rt2_spec) + nrow(rt3_spec))

    while(line_number <= length(raw_data)) {
        current_line <- raw_data[line_number]

        rt1 <- parse_row_with_spec(rt1_spec, current_line)
        
        line_number <- line_number + 1    
        current_line <- raw_data[line_number]

        ## update to read rt2 in loop for number of rt2s
        ## currently not modified yet, this only works for the single rt2 file

        if (is.data.frame(rt2_spec)) {
            
            rt2 <- parse_row_with_spec(rt2_spec, current_line)
            
        } else {
            
            ## update for case where there are multiple rt2 specs
            
        }

        njobs <- rt1['jobs'] %>% parse_integer()


        for (job in 1:njobs) {
            
            line_number <- line_number + 1         
            current_line <- raw_data[line_number]
            rt3 <- parse_row_with_spec(rt3_spec, current_line)
                   
            entries[row_number, ] <- c(rt1, rt2, rt3)
            row_number <- row_number + 1
            
        }
        line_number <- line_number + 1
        
    }
    
    colnames(entries) <- c(pull(rt1_spec, name), pull(rt2_spec, name), pull(rt3_spec, name))
    entries <- entries %>%
        as_tibble() %>%
        filter_all(any_vars( . != ''))

    ## extra modification code for type conversion, this is optional

    if (try_convert) {
        entries %>% 
            mutate(
                dob = parse_date(dob, '%m%d%Y'),
                hire = parse_date(hire, '%m%d%Y'),
                start = parse_date(start, '%m%d%Y'),
                stop = parse_date(stop, '%m%d%Y'),
                jobs = parse_integer(jobs),
                ) %>%
            mutate_if(is.character, parse_guess)
    }
    return(entries)
}

## file <- '/Volumes/Users/Projects/Chloroprene/CommonFiles/OCMAP Files/CD-M_OF5_DEIDENTIFIED.V03'
## entries <- read_ocmap_file(file, rt2_spec = rt2_spec, rt3_spec = rt3_spec)

    


sample_spec <- list(
 ## default rt1 type
rt1_spec = tribble(
    ~name, ~start, ~stop,
    'id', 1, 35,
    'sex', 36, 36,
    'race', 37, 37,
    'cod', 38, 41,
    'dob', 42, 49,
    'hire', 50, 57,
    'start', 58, 65,
    'stop', 66, 73,
    'jobs', 74, 76,
    'unknown_race', 77, 77,
    'vital_status', 78, 78,
    'rate_file_link',  79, 81,
    'race2', 82, 82
),
rt2_spec = tribble(
    ~name, ~start, ~stop,
    'plant_code', 1, 2,
    'latency', 4, 9,
    'duration_of_employment', 11, 16,
    'yob', 18, 21,
    'yoh', 23, 26,
    'yot', 28, 31,
    'yod', 33, 36,
    'aah', 38, 40,
    'raceunk', 42, 42,
    'cd_duration_of_exposure', 44, 49,
    'cd_time_since_1st', 51, 56,
    'cd_cum_exposure', 58, 68,
    'cd_average_intensity', 70, 77,
    'vc_duration_of_exposure', 81, 86,
    'vc_time_since_1st', 88, 93,
    'vc_cum_exposure', 95, 105,
    'vc_average_intensity', 107, 114,
    'cdnovc_duration_of_exposure', 116, 121,
    'cdnovc_time_since_1st', 123, 128,
    'cdnovc_cum_exposure', 130, 135,
    'cdnovc_average_intensity', 137, 14
    ),
rt3_spec <- tribble(
    ~name, ~start, ~stop,
    'job_start_date', 1, 8,
    'job_end_date', 9, 16,
    'job_descriptor', 18, 18,
    'job_duration', 21, 29,
    'job_exposure_cd', 30, 38,
    'job_exposure_vc', 39, 47,
    'job_exposure_nocdnovc', 48, 56,
    'job_exposure_cdnovc', 57, 65,
    'job_exposure_nocdvc', 66, 74,
    'job_exposure_cdvc', 75, 83,
    'job_exposure_cd_presencevc', 84, 92,
    'job_exposure_cd_absencevc', 93, 101,
    'job_duration_cd', 102, 110,
    'job_duration_vc', 111, 119,
    'job_duration_bc', 120, 128,
    'job_duraction_wc', 129, 137
))




library(tidyverse)
library(janitor)
library(readxl)
library(mipfp)
source('table2.R')

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


