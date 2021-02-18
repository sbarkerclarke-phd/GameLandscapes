library(ggmuller)

#Edge generation
edges3 <- data.frame(Parent = paste0("clone_", 
                                     LETTERS[c(1,2,2,3)]), 
                     Identity = paste0("clone_", LETTERS[2:5]))


# a function for generating exponential growth curves:
pop_seq <- function(gens, lambda, start_gen) c(rep(0, start_gen),
                                               exp(lambda * gens[0:(length(gens) - start_gen)]))
lambda <- 0.1 # baseline fitness
gens <- 0:150 # time points
fitnesses <- c(1, 2, 2.2, 2.5, 3, 3.2, 3.5, 3.5, 3.8) # relative fitnesses of genotypes

#Population 
pop3 <- data.frame(Generation = rep(1:150, 5),
                   Identity = paste0("clone_", rep(LETTERS[1:5], each = 150)),
                   Population = c(rep(c(100,90,80),each=50),rep(0,150),rep(0,150),
                                  c(0:29,rep(c(30,35,30,0),each=30)), c(0:29, 29:0, rep(0,30), 0:29, 29:0)))
pop3 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/pop3.csv", colClasses = c("integer", "integer", "character", "double"))
edges3 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/edges3.csv", colClasses = "character")

pop3$X= NULL
edges3$X = NULL
df <- ggmuller::get_Muller_df(edges = edges3, pop_df = pop3)
df <- df[1:200,]
ggmuller::Muller_plot(df, add_legend = TRUE)


pop1 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/pop1.csv", colClasses = c("integer", "integer", "character", "double"))
edges1 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/edges1.csv", colClasses = "character")

pop1$X= NULL
edges1$X = NULL
df <- ggmuller::get_Muller_df(edges = edges1, pop_df = pop1)
df <- df[1:200,]
ggmuller::Muller_plot(df, add_legend = TRUE)


pop2 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/pop2.csv", colClasses = c("integer", "integer", "character", "double"))
edges2 <- read.csv("C:/Users/x_x_s/Documents/ScottLab/GameLandscapes/edges2.csv", colClasses = "character")

pop2$X= NULL
edges2$X = NULL
df <- ggmuller::get_Muller_df(edges = edges2, pop_df = pop2)
df <- df[1:200,]
ggmuller::Muller_plot(df, add_legend = TRUE)
