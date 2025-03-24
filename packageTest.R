sim.data <- sim.marx(c('t',3,0), c('t',2,1), 200, c(0.2,0.3), 0.7, 0.5)
sim.data <- sim.marx(c('t',3,0), NULL, 200, c(0.2,0.3), 0.7, 0)

ts.plot(sim.data$x)

selection.lag(sim.data$y, sim.data$x, 6)
