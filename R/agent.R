library(tidyverse)
library(R6)


## finite time interval state space model


## in agent context
## state(t0) -> state(t0 + interval)

purpose <- 'testing'

evolve <- function(interval){
    params <- self$parameterize(purpose)
    
}



agent <- R6Class('agent',
                 list(
                     evolve = function() {
                         self
                     },
                     receive = function(message){
                         decision = cases(message, agent_message_table, self)
                         decision()
                     }))

cases <- function(message, agent_message_table, self){
    if (message == 'do nothing'){
        function(){}
    }
    else {
        
    }
}
