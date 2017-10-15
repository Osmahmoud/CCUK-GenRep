source("R/Rfuns.R")
# options(MLwiN_path="C:/Program Files/MLwiN v3.01/")
Hubs = c("Cleft Net East" = "CAM", "Leeds" = "LEEDS", "Liverpool" = "LIV", "Manchester" = "MAN", "Midlands" = "MIDS",
         "Northern Ireland" = "N.IRE", "North Thames" = "N.THAMES", "Newcastle" = "NEWC", "South Thames" = "S.THAMES",
         "Scotland" = "SCOT",  "Spires" = "SPIRES", "South West" = "SWSW", "Trent" = "TRENT")

for (i in 1:length(Hubs)){
  tryCatch({
    genReport(Hub = names(Hubs[i]), code.hub = Hubs[i], Overwrite = TRUE)
  }, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
}
