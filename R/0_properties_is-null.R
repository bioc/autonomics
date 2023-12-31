# @rdname is_null
# @export
is_not_null <- function(x, .xname = get_name_in_parent(x))
{
  if(is.null(x))
  {
    return(false("%s is NULL.", .xname))
  }
  TRUE
}

#' Is the input (not) null?
#'
#' Checks to see if the input is (not) null.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_null} wraps \code{is.null}, providing more 
#' information on failure. \code{is_not_null} returns \code{TRUE} in
#' the opposite case.  The \code{assert_*} functions return nothing but
#' throw an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @seealso \code{\link[base]{is.null}}.
#' @examples
#' # Predicate for NULL. 
#' is_null(NULL)
#' is_null(c())
#' 
#' # Atomic vectors of length zero are not NULL!
#' is_null(numeric())
#' # ... and neither is NA
#' is_null(NA)
#' 
#' # The opposite check
#' is_not_null(NULL)
#' is_not_null(c())
#' is_not_null(numeric())
#' 
#' # These checks should pass
#' assert_is_null(NULL)
#' assert_is_null(c())
#' assert_is_not_null(NA)
#' 
#' # This should fail
#' assertive.base::dont_stop(assert_is_null(NaN))
#' @author Richard Cotton
#' @importFrom assertive.base safe_deparse
#' @noRd
is_null <- function(x, .xname = get_name_in_parent(x))
{
  if(!is.null(x))
  {
    return(false("%s is not NULL; its value is %s.", .xname, safe_deparse(x)))
  }
  TRUE
}
