skull = read.csv('C:/Users/daess/desktop/skull.csv')
head(skull)

install.packages('plotly')
library(plotly)

plot_ly(x = ~x1,
        y = ~x2,
        z = ~x3, data = skull) %>%
  add_markers(color = ~group)
head(skull)

#confusion matrix
a = sample(c(0,1), 10, replace = T) #true
b = sample(c(0,1), 10, replace = T) #predict

tb = table(true = a, pred = b)

install.packages('caret')
library(caret) #confusion matrix 함수 포함 패키지지
install.packages('e1071')
library(e1071)
tt = confusionMatrix(tb)
tt$overall
tt$byClass
