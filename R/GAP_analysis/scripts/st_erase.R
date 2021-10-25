#' st_erase
#'
#' @description Removes intersection from two spacial objects. This is useful when removing water areas from TIGER county files
#' @param x sf dataframe
#' @param y sf datafram
#'
#' @return
#' @export
#'
#'
st_erase <- function(x, y) {
  #Check coordinate reference system
  if(st_crs(x) != st_crs(y)){
    y <- st_transform(y, st_crs(x))
  }
  sf::st_difference(x, st_union(y))
}