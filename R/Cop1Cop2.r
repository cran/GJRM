Cop1Cop2 <- function(BivD){

if(BivD == "C0C90"   ){Cop1 <- "C0";   Cop2 <- "C90" }
if(BivD == "C0C270"  ){Cop1 <- "C0";   Cop2 <- "C270"}
if(BivD == "C180C90" ){Cop1 <- "C180"; Cop2 <- "C90" }
if(BivD == "C180C270"){Cop1 <- "C180"; Cop2 <- "C270"}

if(BivD == "GAL0GAL90"   ){Cop1 <- "GAL0";   Cop2 <- "GAL90" }
if(BivD == "GAL0GAL270"  ){Cop1 <- "GAL0";   Cop2 <- "GAL270"}
if(BivD == "GAL180GAL90" ){Cop1 <- "GAL180"; Cop2 <- "GAL90" }
if(BivD == "GAL180GAL270"){Cop1 <- "GAL180"; Cop2 <- "GAL270"}

if(BivD == "J0J90"   ){Cop1 <- "J0";   Cop2 <- "J90" }
if(BivD == "J0J270"  ){Cop1 <- "J0";   Cop2 <- "J270"}
if(BivD == "J180J90" ){Cop1 <- "J180"; Cop2 <- "J90" }
if(BivD == "J180J270"){Cop1 <- "J180"; Cop2 <- "J270"}

if(BivD == "G0G90"   ){Cop1 <- "G0";   Cop2 <- "G90" }
if(BivD == "G0G270"  ){Cop1 <- "G0";   Cop2 <- "G270"}
if(BivD == "G180G90" ){Cop1 <- "G180"; Cop2 <- "G90" }
if(BivD == "G180G270"){Cop1 <- "G180"; Cop2 <- "G270"}

list(Cop1 = Cop1, Cop2 = Cop2)

}


