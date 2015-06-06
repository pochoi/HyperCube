pkgname <- "HyperCube"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('HyperCube')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("bye")
### * bye

flush(stderr()); flush(stdout())

### Name: bye
### Title: Bye world
### Aliases: bye

### ** Examples

## Say Bye!
bye("Taeyen")



cleanEx()
nameEx("hello")
### * hello

flush(stderr()); flush(stdout())

### Name: hello
### Title: Hello world
### Aliases: hello

### ** Examples

## Say hello!
hello("Po")



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
