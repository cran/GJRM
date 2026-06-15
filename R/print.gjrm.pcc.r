print.gjrm.pcc <- function(x, ...) {

    show <- function(obj, label) {
        if (!is.null(obj)) {
            cat("\n---- ", label, " ----", sep = "") # cat("\n--- ", label, " ---\n", sep = "")
            out <- capture.output(print(obj))
            out <- out[out != "[1] NA"]
            cat(out, sep = "\n")
        }
    }

    show(x$b12,   "C12")
    show(x$b13,   "C13")
    show(x$b23,   "C23")
    show(x$b23.1, "C23|1")
    show(x$b13.2, "C13|2")
    show(x$b12.3, "C12|3")

active <- c(
        if (!is.null(x$b12))   "b12",
        if (!is.null(x$b13))   "b13",
        if (!is.null(x$b23))   "b23",
        if (!is.null(x$b23.1)) "b23.1",
        if (!is.null(x$b13.2)) "b13.2",
        if (!is.null(x$b12.3)) "b12.3"
    )

    obj_name <- deparse(substitute(x))

cat(
  "\nFor convergence diagnostics, plots and further analysis,\n use the bivariate component objects: ",
  paste(active, collapse = ", "),
  "\n\n",
  sep = ""
)



    invisible(NULL)
}
