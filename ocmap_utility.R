library(tidyverse)
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
    'dslw', 112, 117,
    'plantev', 119, 122,
    'eagedays', 124, 129,
    'latdays', 131, 136,
    'nosel', 138, 140,
    'noexp', 142, 144,
    'noavg', 146, 148,
    'plant', 150, 165,
    'doe', 166, 181,
    'aah', 198, 213,
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

## calling syntax
data <- read_fwf(file, fwf_positions(start = spec$start, end = spec$stop, col_names = spec$name))


## standard preprocessing
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

