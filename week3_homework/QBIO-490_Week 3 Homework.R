#Week 3 Homework

###Exercise 1.1
#1
attenu[is.na(attenu$station),]

#2
attenu_cleaned = attenu[!is.na(attenu$station),]
attenu_cleaned

#3
head(attenu_cleaned)
dim(attenu_cleaned)

###Exercise 1.2
#1
Theoph_2 = Theoph

#2
median(Theoph_2$Dose) #4.53

#3
Theoph_2$Dose_Class = ifelse(Theoph_2$Dose >= median(Theoph_2$Dose), "high", "low")

#4
head(Theoph_2)
dim(Theoph_2)

###Exercise 1.3
#1
starbucks = read.csv("starbucks.csv")

#2
is.na(starbucks)
is_row_empty = rowSums(is.na(starbucks))
length(is_row_empty) #same length is n rows in starbucks
nrow(starbucks)
starbucks_cleaned = starbucks[is_row_empty == 0,]

#3
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, ylab = "Calories", xlab = "Carbohydrates (g)")

#4
starbucks_cleaned[starbucks_cleaned$Calories == max(starbucks_cleaned$Calories),]$Drink

#5
starbucks_cleaned$is_highest_fat = ifelse(starbucks_cleaned$Fat == max(starbucks_cleaned$Fat), TRUE, FALSE)
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, ylab = "Calories", xlab = "Carbohydrates (g)", 
     col = factor(starbucks_cleaned$is_highest_fat))

#6
library(ggplot2)
ggplot(data=starbucks_cleaned, aes(x = Carb, y = Calories, col = Fat)) + 
  geom_point() +
  scale_colour_gradient(low="green", high="red")

###Exercise 1.4
#1
baseball = read.csv("Batting.csv")

#2
baseball[baseball$HR >= 3,]

#3
plot(baseball$yearID, baseball$HR, xlab = "Year", ylab = "Number of Homeruns")

#4
LAA_players = baseball[baseball$teamID == "LAA",]
plot(LAA_players$yearID, LAA_players$HR, xlab = "Year", ylab = "Number of Homeruns")

#5
ATL_or_PIT_players = baseball[baseball$teamID == "ATL" | baseball$teamID == "PIT",]
plot(ATL_or_PIT_players$yearID, ATL_or_PIT_players$HR, xlab = "Year", ylab = "Number of Homeruns", col = factor(ATL_or_PIT_players$teamID))

###Exercise 1.5
#1
easy_plot <- function(x, y, color_data){
  color_data_median = median(color_data)
  color_levels = ifelse(color_data > color_data_median, "high", "low")
  color_levels = factor(color_levels)
  print(color_data_median)
  plot(x, y, col = color_levels, pch = 20)
  cor.test(x, y)
}
#test of easy_plot, points are colored based on x median
easy_plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, starbucks_cleaned$Carb)

#color by fat content
easy_plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, starbucks_cleaned$Fat)

easy_plot(baseball$H, baseball$RBI, baseball$HR)

###Exercise 2.1
head(iris)
dim(iris)

#The iris dataset contains data on the sepal length and width and petal length and width of 3 different type of irises.
#There are 150 observations and 5 features for each observation.


###Exercise 2.2
#The Sepal.Length, Sepal.Width, Petal.Length, and Petal.Width features are continuous variables (numerics). The Species feature is categorical (characters).

###Exercise 2.3
hist(iris$Sepal.Length)
hist(iris$Sepal.Width)
hist(iris$Petal.Length)
hist(iris$Petal.Width)

###Exercise 2.4
iris_copy = iris
mean_sepal_width = mean(iris$Sepal.Width)
sepal_width_compare = ifelse(iris$Sepal.Width > mean_sepal_width, "wide", "narrow")
iris_copy$Sepal_Width_Factor = sepal_width_compare
boxplot(iris_copy$Sepal.Width ~ iris_copy$Sepal_Width_Factor)

###Exercise 2.5
#The setosa species appears to be the most unique out of the three, while the versicolor and the virginica species are more similar to each other.
iris_subset = iris[,1:4]
pairs(iris_subset, col = factor(iris$Species))


###Exercise 3.1
install.packages("TCGAbiolinks")
library(TCGAbiolinks)
