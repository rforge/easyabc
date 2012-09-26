progressBar <- function (min = 0, max = 1, initial = 0, text = "", char = "=", width = NA, 
    title, label, style = 1, file = "") 
{
    if (!identical(file, "") && !(inherits(file, "connection") && 
        isOpen(file))) 
        stop("'file' must be \"\" or an open connection object")
    if (!style %in% 1L:3L) 
        style <- 1
    .val <- initial
    .killed <- FALSE
    .nb <- 0L
    .pc <- -1L
    nw <- nchar(char, "w")
    if (is.na(width)) {
        width <- getOption("width")
        if (style == 3L) 
            width <- width - 10L
        width <- trunc(width/nw)
    }
    if (max <= min) 
        stop("must have max > min")
    up <- function(value, text) {
        if (!is.finite(value) || value < min || value > max) 
            return()
        .val <<- value
        nb <- round(width * (value - min)/(max - min))
        pc <- round(100 * (value - min)/(max - min))
        if (nb == .nb && pc == .pc) 
            return()
        cat(paste(c("\r  |", rep.int(" ", nw * width + 6)), collapse = ""), file = file)
        cat(paste(c("\r  |", rep.int(char, nb), rep.int(" ", 
            nw * (width - nb)), sprintf("| %3d%% %s", pc,text)), collapse = ""), file = file)
        flush.console()
        .nb <<- nb
        .pc <<- pc
    }
    getVal <- function() .val
    kill <- function() if (!.killed) {
        cat("\n", file = file)
        flush.console()
        .killed <<- TRUE
    }
    up(initial, text)
    structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}

updateProgressBar <- function (pb, value, text = "") 
{
    if (!inherits(pb, "txtProgressBar")) 
        stop("'pb' is not from class \"txtProgressBar\"")
    oldval <- pb$getVal()
    pb$up(value,text)
    invisible(oldval)
}

# Example
pb <- progressBar(width=30)
simus = 2*c(0.9, 1.2, 0.8, 0.9, 1);
duration = 0;
for(i in 1:length(simus)) {
  start = Sys.time()
  Sys.sleep(simus[i]);
  duration = duration + difftime(Sys.time(), start, unit="secs")
  text="";
  if (i==length(simus)) {
    text = paste("Completed  in",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"                                              ");
  } else {
    text = paste("− Time elapsed:",format(.POSIXct(duration, tz="GMT"), "%H:%M:%S"),"− Estimated time remaining:",format(.POSIXct(duration/i*(length(simus)-i), tz="GMT"), "%H:%M:%S"));
  }
  updateProgressBar(pb, i/length(simus), text)
}
close(pb)

