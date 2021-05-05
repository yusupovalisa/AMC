setwd('/home/npav/new_code/ijf2020/')

require(data.table)
require(ggplot2)
require(dplyr)

# ABRUPT CHANGE
X <- fread("afSims/AFDLM_sim_abrupt.csv")
# Plot evolution of estimated coefficients
pp <- ggplot(X) + 
#  scale_color_manual(values=c("red", "green","black","blue","yellow")) +
#  scale_fill_manual(values=c("red", "green", "black","blue","yellow")) +
  geom_ribbon(aes(x=time, ymin=q1th1, ymax=q3th1, fill="red"), alpha=0.15, show.legend =FALSE)+
  geom_line(aes(x=time, y=medth1,color="red"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  geom_line(aes(x=time, y=th1,color="red"), size=0.5, linetype = "solid", show.legend = FALSE) +
  #
  geom_ribbon(aes(x=time, ymin=q1th2, ymax=q3th2, fill="green"), alpha=0.15, show.legend =FALSE)+
  geom_line(aes(x=time, y=medth2,color="green"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  geom_line(aes(x=time, y=th2,color="green"), size=0.5, linetype = "solid", show.legend = FALSE) +
  #
  geom_ribbon(aes(x=time, ymin=q1th3, ymax=q3th3, fill="black"), alpha=0.15, show.legend =FALSE)+
  geom_line(aes(x=time, y=medth3,color="black"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  geom_line(aes(x=time, y=th3,color="black"), size=0.5, linetype = "solid", show.legend = FALSE) +
  #
  geom_ribbon(aes(x=time, ymin=q1th4, ymax=q3th4, fill="blue"), alpha=0.15, show.legend =FALSE)+
  geom_line(aes(x=time, y=medth4,color="blue"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  geom_line(aes(x=time, y=th4,color="blue"), size=0.5, linetype = "solid", show.legend = FALSE) +
  #
  geom_ribbon(aes(x=time, ymin=q1th5, ymax=q3th5, fill="yellow"), alpha=0.15, show.legend =FALSE)+
  geom_line(aes(x=time, y=medth5,color="yellow"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  geom_line(aes(x=time, y=th5,color="yellow"), size=0.5, linetype = "solid", show.legend = FALSE) +
  #
  theme_bw() + 
  theme(text = element_text(size=12)) + theme(plot.background = element_blank(), panel.grid.major = element_blank() ) +
  theme(panel.grid.minor = element_blank() ) +
  ylab(expression("True and Estimated Parameters" ~ ~ (theta[t]))) + xlab("Time") +
  scale_y_continuous(breaks=c(seq(-2,3,1)), limits=c(-2.5,3.5)) +
  scale_x_continuous(breaks=c(seq(0,1000,100)), limits=c(0,1000)) 
  
print(pp)
ggsave(filename="abrupt_coeff.pdf", plot=pp, width=10, height=8, units="in")

# Plot evolution of forgetting factor
pp <- ggplot(X) + 
  geom_line(aes(x=time, y=medL), size=1, linetype="solid", color="blue",show.legend = FALSE) +
  geom_ribbon(aes(x=time, ymin=q1L, ymax=q3L), fill="blue", alpha=0.15, show.legend =FALSE) +
  theme_bw() + 
  theme(text = element_text(size=12)) + theme(plot.background = element_blank()) +#, panel.grid.major = element_blank() ) +
  theme(panel.grid.minor = element_blank() ) +
  xlab("Time") +  ylab(expression("Forgeting Factor" ~ (lambda))) +
  scale_x_continuous(breaks=c(seq(0,1000,100)), limits=c(0,1000)) 
print(pp)
ggsave(filename="abrupt_lambda.pdf", plot=pp, width=10, height=8, units="in")



###### GRADUAL CHANGE
X <- fread("afSims/AFDLM_sim_gradual.csv")
pp <- ggplot() + 
  geom_line(data=X[lambda==0.99,], aes(x=time, y=medL), size=1, linetype="solid", color="blue",show.legend = FALSE) +
  geom_ribbon(data=filter(X,lambda==0.99), aes(x=time, ymin=q1L, ymax=q3L), fill="blue", alpha=0.15, show.legend =FALSE) +
  geom_hline(yintercept = 0.99, color="blue", linetype="dashed")  +
  #
  geom_line(data=filter(X,lambda==0.97), aes(x=time, y=medL), size=1, linetype="solid", color="green",show.legend = FALSE) +
  geom_ribbon(data=filter(X,lambda==0.97), aes(x=time, ymin=q1L, ymax=q3L), fill="green", alpha=0.15, show.legend =FALSE) +
  geom_hline(yintercept = 0.97, color="green", linetype="dashed")  +
  #
  geom_line(data=filter(X,lambda==0.95), aes(x=time, y=medL), size=1, linetype="solid", color="red",show.legend = FALSE) +
  geom_ribbon(data=filter(X,lambda==0.95), aes(x=time, ymin=q1L, ymax=q3L), fill="red", alpha=0.15, show.legend =FALSE) +
  geom_hline(yintercept = 0.95, color="red", linetype="dashed")  +
  theme_bw() + 
  theme(text = element_text(size=12)) + theme(plot.background = element_blank()) +#, panel.grid.major = element_blank() ) +
  theme(panel.grid.minor = element_blank() ) +
  ylab(expression("Forgeting Factor" ~ (lambda))) + xlab("Time") +
  scale_x_continuous(breaks=c(seq(0,1000,100)), limits=c(0,1000)) 
print(pp)
ggsave(filename="gradual_lambda.pdf", plot=pp, width=10, height=8, units="in")



# STATIC 
X <- fread("afSims/AFDLM_sim_static.csv")
# Plot evolution of estimated coefficients
pp <- ggplot(X[time<=30,]) + 
  geom_ribbon(aes(x=time, ymin=q1th1, ymax=q3th1, fill="red"), alpha=0.15, show.legend =FALSE)+
  geom_line(aes(x=time, y=medth1,color="red"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  geom_line(aes(x=time, y=th1,color="red"), size=0.5, linetype = "solid", show.legend = FALSE) +
  #
  geom_ribbon(aes(x=time, ymin=q1th2, ymax=q3th2, fill="green"), alpha=0.15, show.legend =FALSE)+
  geom_line(aes(x=time, y=medth2,color="green"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  geom_line(aes(x=time, y=th2,color="green"), size=0.5, linetype = "solid", show.legend = FALSE) +
  #
  geom_ribbon(aes(x=time, ymin=q1th3, ymax=q3th3, fill="black"), alpha=0.15, show.legend =FALSE)+
  geom_line(aes(x=time, y=medth3,color="black"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  geom_line(aes(x=time, y=th3,color="black"), size=0.5, linetype = "solid", show.legend = FALSE) +
  #
  geom_ribbon(aes(x=time, ymin=q1th4, ymax=q3th4, fill="blue"), alpha=0.15, show.legend =FALSE)+
  geom_line(aes(x=time, y=medth4,color="blue"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  geom_line(aes(x=time, y=th4,color="blue"), size=0.5, linetype = "solid", show.legend = FALSE) +
  #
  geom_ribbon(aes(x=time, ymin=q1th5, ymax=q3th5, fill="yellow"), alpha=0.15, show.legend =FALSE) +
  geom_line(aes(x=time, y=medth5,color="yellow"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  geom_line(aes(x=time, y=th5, color="yellow"), size=0.5, linetype = "solid", show.legend = FALSE) +
  #
  theme_bw() + 
  theme(text = element_text(size=12)) + theme(plot.background = element_blank(), panel.grid.major = element_blank() ) +
  theme(panel.grid.minor = element_blank() ) +
  ylab(expression("True and Estimated Parameters" ~ ~ (theta))) + xlab("Time") +
  scale_y_continuous(breaks=c(seq(-2,3,1)), limits=c(-3,3.5)) 
  #scale_x_continuous(breaks=c(seq(0,1000,100)), limits=c(0,1000)) 

print(pp)
ggsave(filename="static_coeff.pdf", plot=pp, width=10, height=8, units="in")

# Plot evolution of forgetting factor
pp <- ggplot(X) + 
  geom_line(aes(x=time, y=medL), size=1, linetype="solid", color="blue",show.legend = FALSE) +
  geom_ribbon(aes(x=time, ymin=q1L, ymax=q3L), fill="blue", alpha=0.15, show.legend =FALSE) +
  theme_bw() + 
  theme(text = element_text(size=12)) + theme(plot.background = element_blank()) +#, panel.grid.major = element_blank() ) +
  theme(panel.grid.minor = element_blank() ) +
  ylab(expression("Forgeting Factor" ~ (lambda))) +
  scale_x_continuous(breaks=c(seq(0,1000,100)), limits=c(0,1000))  +
  scale_y_continuous(limits=c(0.95,1)) 
print(pp)
ggsave(filename="static_lambda.pdf", plot=pp, width=10, height=8, units="in")


##################### MODEL COMBINATION EXPERIMENTS
# Persistence
X <- fread("mcSims/ConfHedge_persistence.csv")
X <- X[,c("phi","time","medW1","medW2","medW3")]
names(X) <- c("phi","time","w1","w2","w3")
X <- data.table::melt(X, id.vars=c("time","phi"))

hp <- ggplot(data = X) +
  geom_tile(aes(x=time,y=phi, fill=value)) +
  scale_fill_gradient(low="white", high="black", guide = "colourbar", breaks=seq(0,1,length.out = 11)) +
  scale_y_continuous(breaks=c(0,1,seq(0,1,length.out = 11)), limits = c(0.02,0.99), name =  expression(phi)) +
  scale_x_continuous(breaks=c(50,300,seq(0,300,by=50)), name="time") +
  theme_bw() +
  theme(text = element_text(size=22), legend.text = element_text(size=22),
        plot.background = element_blank(),
        panel.grid.major =  element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_rect(size=0),
        legend.title = element_text(color="white"),
        legend.key.height = unit(2.5, "cm")
  ) +
  facet_wrap(variable ~ ., ncol=3)
print(hp)
ggsave(filename="mcSims/persistence.pdf", plot=hp, width=25, height=8, units="in")

# Abrupt Change
X <- fread("mcSims/abrupt_change.csv")
names(X)[1] <- "algorithm"
Y <- data.table::melt(X, id.vars=c("time","algorithm"), measure.vars=c("w1","w2","w3"), variable.name="weight", value.name="median")
Y1 <- data.table::melt(X, id.vars=c("time","algorithm"), measure.vars=c("q1w1","q1w2","q1w3"), variable.name="weight", value.name="q1")
Y3 <- data.table::melt(X, id.vars=c("time","algorithm"), measure.vars=c("q3w1","q3w2","q3w3"), variable.name="weight", value.name="q3")
Y <- cbind(Y, Y1[,"q1"], Y3[,"q3"])
Y$algorithm <- factor(Y$algorithm, levels = c("BMA","DMA99","DMA95","OPP","ConfHedge"))
pp <- ggplot(Y) + 
  geom_ribbon(aes(x=time, ymin=q1, ymax=q3, fill=algorithm), alpha=0.15, show.legend = FALSE)+
  geom_line(aes(x=time, y=median,color=algorithm, linetype=algorithm), size=0.5, show.legend = TRUE) +
  #
  # geom_ribbon(aes(x=time, ymin=q1th2, ymax=q3th2, fill="green"), alpha=0.15, show.legend =FALSE)+
  # geom_line(aes(x=time, y=th2,color="green"), size=0.5, linetype = "solid", show.legend = FALSE) +
  # #
  # geom_ribbon(aes(x=time, ymin=q1th3, ymax=q3th3, fill="black"), alpha=0.15, show.legend =FALSE)+
  # geom_line(aes(x=time, y=medth3,color="black"), size=0.5, linetype = "dotted", show.legend =FALSE) +
  # geom_line(aes(x=time, y=th3,color="black"), size=0.5, linetype = "solid", show.legend = FALSE) +
  # #
  # geom_ribbon(aes(x=time, ymin=q1th4, ymax=q3th4, fill="blue"), alpha=0.15, show.legend =FALSE)+
  # geom_line(aes(x=time, y=th4,color="blue"), size=0.5, linetype = "solid", show.legend = FALSE) +
  # #
  # geom_ribbon(aes(x=time, ymin=q1th5, ymax=q3th5, fill="yellow"), alpha=0.15, show.legend =FALSE)+
  # geom_line(aes(x=time, y=th5,color="yellow"), size=0.5, linetype = "solid", show.legend = FALSE) +
  # #
  geom_vline(xintercept = c(100,200), linetype="dashed") +
  theme_bw() + xlab("time") + ylab("") +
  theme(text = element_text(size=18), legend.position = "bottom",
        legend.text = element_text(size=18), legend.key.width = unit(1.,"cm")) +
  #theme(text = element_text(size=12)) + theme(plot.background = element_blank(), panel.grid.major = element_blank() ) +
  #theme(panel.grid.minor = element_blank() ) +
  facet_grid(weight ~ .)
  # ylab(expression("True and Estimated Parameters" ~ ~ (theta[t]))) + xlab("Time") +
  # scale_y_continuous(breaks=c(seq(-2,3,1)), limits=c(-2.5,3.5)) +
  # scale_x_continuous(breaks=c(seq(0,1000,100)), limits=c(0,1000)) 
print(pp)
ggsave(filename="mcSims/persistence.pdf", plot=pp, width=10, height=8, units="in")

pp <- ggplot(X) + geom_line(aes(x=time,y=))