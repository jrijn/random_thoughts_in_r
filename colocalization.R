# This function can be used to check e.g. colocalization of object locations from ImageJ>Analyze particles

require(tidyverse)
require(pdist)

colocalization <- function(.data,               # The input data.frame, supports %>% pipe.
                           id,                  # A column used as unique group identifier, e.g. one field-of-view.
                           var,                 # A column identifying two object groups to check for colocalization. E.g. channel.
                           check,               # A string of the reference objects. E.g. the reference channel 'mcherry'.
                           against,             # A string of the objects checked against. E.g. to check for 'AF488' colocalization.
                           ...,                 # Unquoted column names with the dimensions. E.g. X, Y(, Z, ...).
                           coloc.dist = 2){     # The minimum distance to be considered colocalized. 
  
  dims <- tidyselect::eval_select(expr(c(...)), .data)
  ids <- .data %>% distinct({{id}}) %>% pull({{id}})
  out <- data.frame()
  for (i in ids) {
    df.i <- .data %>% filter({{id}} == i)
    df.check <- df.i %>% filter({{var}} == check) %>% select(all_of(dims))
    df.against <- df.i %>% filter({{var}} == against) %>% select(all_of(dims))
    df.pdist <- as.matrix(pdist(df.check, df.against))
    df.out <- df.i %>%
      filter({{var}} == check) %>%
      mutate(rownumber = 1:n()) %>%
      rowwise() %>%
      mutate(coloc = any(df.pdist[rownumber,] < coloc.dist)) %>%
      select(-rownumber)
    out <- bind_rows(out, df.out)
  }
  return(out)
}

test <- data.frame(
  condition = rep(c("A", "B"), times = 3, each = 10),
  replicate = rep(1:3, times = 1, each = 20),
  channel = rep(c("channel1", "channel2"), times = 6, each = 5),
  X = runif(60, min = 0, max = 100),
  Y = runif(60, min = 0, max = 100)
)
test

out <- test %>%
  unite(ID, condition, replicate, remove = FALSE) %>%
  colocalization(ID, channel, "channel1", "channel2", X, Y, coloc.dist = 20)
out
