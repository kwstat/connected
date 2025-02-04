
The `apple_production` data is taken from USApple Industry Outlook reports for 2021-2023.
https://usapple.org/wp-content/uploads/2021/08/USAppleIndustryOutlook2021.pdf
https://usapple.org/wp-content/uploads/2022/08/USAPPLE-INDUSTRYOUTLOOK-2022.pdf
https://usaa.memberclicks.net/assets/2023/Outlook2023/USAPPLE-OutlookReport-2023.pdf

The values are the *estimated* number of million bushels of apples produced, broken out by year, state, and variety.

apples <- read.csv("c:/drop/rpack/connected/data-raw/apple_production.csv")
head(apples)
library(reshape2)
apples <- melt(apples, id.var=c("year","variety"))
apples <- rename(apples, region=variable, bushels=value)
apples <- mutate(apples, bushels=bushels/1000000)
head(apples)
apples %>% filter(!(region %in% c("OTHER","US"))) %>%
  mutate(lbushels=log(bushels)) %>% 
  con_view(lbushels ~ variety * year:region)
