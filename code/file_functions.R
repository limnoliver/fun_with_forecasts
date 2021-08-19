# functions to make it happen
return_da_files <- function(file_location, model) {
  all_files <- list.files(file_location, full.names = TRUE)
  mod_files <- grep('\\.nc$', all_files, value = TRUE)
  mod_strings <- tibble(model_strings = c('_no_DA_', '\\]_DA_', 'persistence'),
                        model_name = c('no_DA', 'DA', 'persistence')) %>%
    filter(model_name %in% model)
  mod_files_out <- grep(paste(mod_strings$model_strings, collapse = "|"), mod_files, value = TRUE)
  
  return(mod_files_out)
}

extract_and_bind <- function(in_files){
  out_conservation <- map2_dfr(
    .x = in_files[grepl('\\[0_cfs\\]', in_files)],
    .y = 'max_temp',
    .f = ~ nc_forecast_get(.x, var_name = .y, issue_times = 'max')) %>%
    filter(!is.na(max_temp))
  
  return(out_conservation)
}
#' function for returning forecast variables; returns tibble
#'
#' @param nc_file file path to the netcdf file
#' @param var_name forecast variable name you want returned from the netcdf file
#' @param times model times you want returned
#' @param seg_id_nats seg_id_nat segments you want returned
#'
nc_forecast_get = function(nc_file,
                           var_name,
                           issue_times = NULL,
                           times = NULL,
                           seg_id_nats = NULL,
                           ens = NULL){
  
  ncout = nc_open(nc_file)
  
  # output dims [issue_time_dim, time_dim, loc_dim, ens_dim]; position of dimensions following:
  issue_time_pos = 1
  time_pos = 2
  loc_pos = 3
  ens_pos = 4
  
  cur_var = ncout$var[[var_name]]
  varsize = cur_var$varsize
  all_issue_times = c(ncdf_times(nc = ncout, timename = 'issue_time'))
  max_issue_time = max(all_issue_times)
  all_valid_times = cur_var$dim[[time_pos]]$vals
  all_seg_id_nats = as.character(cur_var$dim[[loc_pos]]$vals)
  all_ens = as.integer(cur_var$dim[[ens_pos]]$vals)
  
  n_dims = cur_var$ndims
  
  # return all values, and then filter
  all_out = array(ncvar_get(nc = ncout, varid = var_name), dim = varsize) %>%
    reshape2::melt(varnames = c('issue_time', 'time', 'seg_id_nat', 'ensemble')) %>%
    mutate(issue_time = all_issue_times[issue_time],
           time = issue_time + as.difftime(all_valid_times[time], units = 'days'),
           seg_id_nat = all_seg_id_nats[seg_id_nat],
           ensemble = all_ens[ensemble]) %>%
    rename(!!var_name := value) %>%
    as_tibble()
  
  if(!is.null(issue_times) & !('max' %in% issue_times)){
    cur_issue_times = as.Date(issue_times)
  }else if ('max' %in% issue_times){
    cur_issue_times = as.Date(max_issue_time)
  } else {
    cur_issue_times = as.Date(all_issue_times)} # return all issue times if NULL
  if(!is.null(times)){
    cur_times = as.integer(times)
  }else{cur_times = as.integer(all_valid_times)}
  if(!is.null(seg_id_nats)){
    cur_seg_id_nats = as.character(seg_id_nats)
  }else{cur_seg_id_nats = as.character(all_seg_id_nats)}
  if(!is.null(ens)){
    cur_ens = as.integer(ens)
  }else{cur_ens = as.integer(all_ens)}
  out = dplyr::filter(all_out,
                      issue_time %in% cur_issue_times,
                      seg_id_nat %in% cur_seg_id_nats,
                      ensemble %in% cur_ens) %>%
    mutate(filter_valid_time = as.integer(time - issue_time)) %>%
    dplyr::filter(filter_valid_time %in% cur_times) %>%
    dplyr::select(-filter_valid_time) %>%
    mutate(model_name = gsub('(.*\\]_)(.*)(_forecast.*)', '\\2', nc_file, perl = TRUE))
  
  nc_close(ncout)
  
  return(out)
}

#' Get time attribute from NetCDF file and convert it to years (or Rdate)
#'
#' This function reads the NetCDF time attribute and converts it to years
#' with decimal fraction for monthly data. For sub-monthly data the conversion
#' is to the Rdate format (only works with standard calendar!).
#'
#' original code from https://rdrr.io/github/jonasbhend/geoutils/src/R/ncdf_times.R
#'
#' @param nc object with link to NetCDF file (from \code{\link{open.ncdf}})
#' @param timename dimension name with time (optional)
#' @param as.Rdate logical, should output be converted to Rdate?
#' @param force logical, force Rdate conversion for monthly times?
#' @param tz time zone for time zone support (experimental)
#'
ncdf_times <- function(nc, timename=NULL, as.Rdate=TRUE, force=TRUE, tz="UTC") {
  ## this function converts netcdf times to the
  ## R date-time format or to the udunits dates
  ## you can choose to switch to uduints format
  ## for gregorian calendars by setting as.Rdate
  ## to FALSE
  ## non-gregorian calendar dates are output using
  ## udunits date format
  
  ## you can force to get back an R date format, even
  ## if the calendar used is not gregorian using
  ## force=T (may return udunits dates if conversion
  ## is not successful)
  
  ## WARNING: time zones are not fully supported yet
  
  ## check whether udunits is available
  .udunitsInclude     <- FALSE
  if (any(.packages() == "udunits") & class(try(utInit(), silent=T)) != "try-error"){
    .udunitsInclude <- TRUE
  }
  if (is.null(timename)){
    timei   <- which(names(nc$dim) %in% c("time", "TIME", "tim", "TIM"))
  } else {
    timei <- which(names(nc$dim) == timename)
  }
  units   <- nc$dim[[timei]]$units
  refdate <- strsplit(units, " ")[[1]]
  refdate <- refdate[grep('-', refdate)]
  ## debug reference date
  if (substr(refdate, nchar(refdate) - 1, nchar(refdate)) == '00') {
    rtmp <- strsplit(refdate, '-')[[1]]
    refdate <- paste(c(rtmp[-length(rtmp)], '01'), collapse='-')
    rm(rtmp)
  }
  vals    <- nc$dim[[timei]]$vals
  tmp     <- ncatt_get(nc, names(nc$dim)[timei], "calendar")
  if (tmp$hasatt) {
    calendar <- tmp$value
  } else {
    calendar <- "standard"
    ## print(paste("Warning: Calendar is missing in", nc$filename))
  }
  if (calendar == "proleptic_gregorian" || calendar == "gregorian") calendar <- "standard"
  
  if (as.Rdate){
    
    if (charmatch("hours since", units, nomatch=0) |
        charmatch("minutes since", units, nomatch=0) |
        charmatch("seconds since", units, nomatch=0) & calendar == 'standard') {
      
      mul <- 1
      ref.txt     <- substr(units, 15,33)
      if (charmatch("minutes", units, nomatch=0)) mul     <- 60
      if (charmatch("hours", units, nomatch=0)) {
        mul     <- 3600
        ref.txt <- substr(units, 13,31)
      }
      
      times       <- vals * mul
      if (nchar(ref.txt) == 19){
        ref         <- as.POSIXct(ref.txt, tz)
      } else {
        ref         <- as.POSIXct(paste(ref.txt, "00", sep=":"), tz)
      }
      time        <- as.Date(ref + times)
      
    } else if (charmatch("days since", units, nomatch=0) & calendar == 'standard'){
      time        <- as.Date(refdate, "%Y-%m-%d") + vals
    } else if (charmatch("days since", units, nomatch=0) &
               calendar %in% c('365_day', 'noleap', '360_day')) {
      if (calendar == '365_day' || calendar == 'noleap'){
        vals <- round(vals/365*365.24*2)/2
        time <- as.Date(refdate, "%Y-%m-%d") + vals
      } else if (calendar == '360_day'){
        vals <- round(vals/360*365.24*2)/2
        time <- as.Date(refdate, "%Y-%m-%d") + vals
      }
    } else if (charmatch("months since", units, nomatch=0)) {
      ref.yr      <- as.numeric(format(as.Date(refdate), '%Y'))
      ref.month   <- as.numeric(format(as.Date(refdate), '%m'))
      ref.day     <- as.numeric(format(as.Date(refdate), '%d'))
      if (is.null(ref.day) | ref.day == 0) ref.day <- 1
      
      month       <- floor((vals+ref.yr*12 + ref.month-1) %% 12) + 1
      year        <- floor((vals+ref.yr*12 + ref.month-1)/12)
      
      time        <- as.Date(ISOdate(year, month, ref.day))
    } else if (charmatch("years since", units, nomatch=0)) {
      unit.tmp    <- paste(strsplit(units, " ")[[1]][3:4], collapse=" ")
      ## ref.yr      <- substr(units, 13,16)
      ## ref.month   <- as.numeric(substr(units, 18,19))
      ref.yr      <- as.numeric(format(as.Date(unit.tmp), "%Y"))
      ref.month   <- as.numeric(format(as.Date(unit.tmp), "%m"))
      if (is.null(ref.month)) ref.month <- 1
      ##ref.day     <- as.numeric(substr(units, 21,22))
      ref.day     <- as.numeric(format(as.Date(unit.tmp), "%d"))
      if (is.null(ref.day)) ref.day <- 1
      
      year        <- floor(vals)
      month       <- floor((vals*12)%%12)
      
      time        <- as.Date(ISOdate(ref.yr + year, ref.month + month, ref.day))
    }  else if (charmatch("day as", units, nomatch=0)) {
      date        <- floor(vals)
      day         <- as.numeric(substr(date, nchar(date)-1, nchar(date)))
      if (all(day > 28)) date <- as.character(as.numeric(date) - max(day, na.rm=T) + 28)
      date        <- paste("000",date, sep="")
      date        <- substr(date, nchar(date)-7, nchar(date))
      time        <- as.Date(date, "%Y%m%d")
    } else {
      stop(paste("Can't deal with calendar", calendar))
    }
    
  } else {
    if (.udunitsInclude){
      time <- utCalendar(vals, units, calendar=calendar, style="array")
      if (force){
        tmp  <- try(ISOdatetime(time$year, time$month, time$day, time$hour,
                                time$minute, time$second, tz), silent=T)
        if (class(tmp)[1] != "try-error") time <- tmp
      }
    } else {
      warning("Package udunits cannot be loaded or initialized via utInit()")
    }
  }
  
  
  return(time)
  
}