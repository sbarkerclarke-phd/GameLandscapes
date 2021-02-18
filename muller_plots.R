library(ggmuller)

#Population 

pop0 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/pop0.csv", colClasses = c("integer", "integer", "character", "double"))
edges0 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/edges0.csv", colClasses = "character")


pop1 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/pop1.csv", colClasses = c("integer", "integer", "character", "double"))
edges1 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/edges1.csv", colClasses = "character")

pop2 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/pop2.csv", colClasses = c("integer", "integer", "character", "double"))
edges2 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/edges2.csv", colClasses = "character")

pop3 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/pop3.csv", colClasses = c("integer", "integer", "character", "double"))
edges3 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/edges3.csv", colClasses = "character")

pop4 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/pop4.csv", colClasses = c("integer", "integer", "character", "double"))
edges4 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/edges4.csv", colClasses = "character")

pop0$X= NULL
edges0$X = NULL
df <- ggmuller::get_Muller_df(edges = edges0, pop_df = pop0)
ggmuller::Muller_plot(df, add_legend = TRUE)

pop1$X= NULL
edges1$X = NULL
df <- ggmuller::get_Muller_df(edges = edges1, pop_df = pop1)
ggmuller::Muller_plot(df, add_legend = TRUE)

pop2$X= NULL
edges2$X = NULL
df <- ggmuller::get_Muller_df(edges = edges2, pop_df = pop2)
ggmuller::Muller_plot(df, add_legend = TRUE)

pop3$X= NULL
edges3$X = NULL
df <- ggmuller::get_Muller_df(edges = edges3, pop_df = pop3)
ggmuller::Muller_plot(df, add_legend = TRUE)


pop4$X= NULL
edges4$X = NULL
df <- ggmuller::get_Muller_df(edges = edges4, pop_df = pop4)
df <- df
ggmuller::Muller_plot(df, add_legend = TRUE)
