
library(ChemmineR)


v2kLines = read.SDFstr("~/runs/v3000/DrugLike-0_2-3K3K_1.v2k.sdf")
v3kLines = read.SDFstr("~/runs/v3000/DrugLike-0_2-3K3K_1.sdf")


message("v3000: ")
system.time(read.SDFset(v3kLines))

message("header time: ", ChemmineR:::v3kTimes$header)
message("atom time: ", ChemmineR:::v3kTimes$atom)
message("atom core time: ", ChemmineR:::v3kTimes$atomCore)
message("bond time: ", ChemmineR:::v3kTimes$bond)
message("data time: ", ChemmineR:::v3kTimes$data)

message("v2000: ")
system.time(read.SDFset(v2kLines))
message("header time: ", ChemmineR:::v2kTimes$header)
message("atom time: ", ChemmineR:::v2kTimes$atom)
message("bond time: ", ChemmineR:::v2kTimes$bond)
message("data time: ", ChemmineR:::v2kTimes$data)


