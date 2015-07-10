# convert excel and other junk to Rds (or csv)

library(haven)
library(readr)

h <- read_csv("../data/regional_data.csv")
glimpse(h)

# Переименуем для удобства:
## ------------------------------------------------------------------------
h <- dplyr::rename(h, y_star = Y, Wy_star = WbY, y_0 = `ln(gdppercappp)`, region = `[EMPTY]`) %>% 
  dplyr::select(-number)

# матрица близости по длине общей границы
W_border <- read_csv("../data/Wb.csv",col_names = FALSE)
W_border <- as.matrix(W_border)
str(W_border)

# матрица близости по расстояниям по автодороге в минус первой
W_road <- read_dta("../data/Russia_Inverted_distance.dta")
W_road <- as.matrix(W_road)
str(W_road)

# нормируем матрицы (сумма по каждой строке должна равняться 1)
apply(W_border, 1, sum)
apply(W_road, 1, sum)

W_road <- W_road/apply(W_road, 1, sum)

nrow(h)
nrow(W_road)
nrow(W_border)

# первые 52 региона западные (в том числе Калининград), а остальные 23 восточные.
h$location <- c(rep("west",52),rep("east",23))

# ищем Калининградскую область
n_kalin <- which(h$region=="Kaliningrad region")
sum(W_border[n_kalin,])

# Она ни с кем не граничит, поэтому при пограничной $W$ 
# нарушается свойство $W\vec{1}=\vec{1}$ и определитель $det(I_{n\times n }-\rho W)$ оказывается отрицательным. 
# А он фигурирует в плотностях. 

h_no_kalin <- h[-n_kalin,]
W_border <- W_border[-n_kalin,-n_kalin]

# сохраняем всё нужное чистое и готовое
write_csv(h_no_kalin,"../data/h_nk.csv")
write_csv(data.frame(W_border),"../data/W_border.csv")

write_csv(h,"../data/h.csv")
write_csv(data.frame(W_road),"../data/W_road.csv")



