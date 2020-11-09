dyn.load("strict_split.so")

str1 <- c("hello,there,", ",kdkd,,", "kdk,kd , akd , lsls,")

for(i in 1:1000000){
    a <- .Call("strict_split", str1, ",")
}

.Call("strict_split", str1, "")

