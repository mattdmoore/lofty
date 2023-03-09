#
#
#
#Author: James Lofty
#script characterises the trajectory and collision charateristcs of spherical microplastic particles based on coordinates (x,y) generated from a high speed camera extracted from the example video
#
#
#

#data inout#####################################################################
#read text file
data = read.csv('data/track 1.txt',  
                sep =' ',                                                       # specify tab-separate                                       
                col.names = c('x', 'y'),                                        # name columns                                        
                na.strings = 'null')                                            #strip NA rows

#raw data processing############################################################
data = data[!is.na(data$x),]                                                    #strip NA rows                                                    
data$y = -data$y + max(data$y)                                                  # flip y and re-center
data = ((data*0.264)*0.80048157)                                                #camera correction (pixels to mm. camera distance correction)

#Variables######################################################################
fs = 90                                                                         #frame rate (s)
MP_d = 5                                                                        #MP diameter (mm)
ks = 1.86                                                                       #roughness height (mm)
susthres = 0.07                                                                 #suspension mode threshold
rollthreshp = ks/MP_d                                                           #rolling mode threshold

#data processing################################################################
#peaks and troughs function
remove_consecutive = function(troughs, peaks)
{
  zipped = rbind(troughs, peaks)
  bad_idx = which.max(diff(c(zipped)) < 0) + 1
  bad_val = zipped[bad_idx]
  if(zipped[2,1] == zipped[2, length(zipped[2,])])
  {
    return(list(
      troughs = troughs, 
      peaks = peaks))
  }
  peaks = peaks[!peaks %in% bad_val]
  troughs = troughs[!troughs %in% bad_val]
  
  # Recurse
  remove_consecutive(troughs, peaks)
}

# Find direction changes
direction = diff(diff(data$y) > 0)                                              

# Find all peaks and troughs
peaks = which(direction == -1) + 1                                              
troughs = which(direction == 1) + 1                                             

# Find the erroneous peaks and troughs
dips = which(diff(direction) == -2) +1
bumps = which(diff(direction) == 2) +1

# Remove erroneous peaks and troughs
peaks = peaks[!peaks %in% bumps]
troughs = troughs[!troughs %in% dips]
peaks = peaks[peaks > min(troughs) & peaks < max(troughs)]

# Remove consecutive peaks and troughs
cleaned = remove_consecutive(troughs, peaks)
troughs = cleaned$troughs
peaks = cleaned$peaks

# angles function
in_out = function(data, trough)
{
  vectors = apply(data[(trough-1):(trough+1),], 2, diff)
  return(apply(vectors, 1, function(x) atan(x[2]/x[1])*(180/pi)))               #inwards and outwards angles equation relartive to flat bed.
}


reld = (data$y[troughs[-1]] - data$y[peaks])/MP_d                               #relative particle drop length 
Z0 = (data$y[peaks]-mean(data$y[troughs]))/MP_d                                 #relative height of grain bed

#raw_results####################################################################
saltations = seq_along(troughs[-1])                                             #saltation number. saltation 1 is from first trough to second
t = diff(troughs)/fs                                                            # time between troughs in seconds
Lp = (diff(data$x[troughs]))                                                    #saltation distance on X axis
Hp = (data$y[peaks] - data$y[troughs[-length(troughs)]])                        #height of saltation from bottom of trough to top of peak
Up = ((diff(data$x[troughs]))/(diff(troughs)/fs)/1000)                          #mean saltation velocity in m/s
angles = sapply(troughs, function(x) in_out(data, x))                           #inwards and outwards collision angle function
Uin = ((data$x[troughs]-data$x[troughs-1])/(1/fs)/1000)                         #inwards streamwise collision velocity
Uout = ((data$x[troughs+1]-data$x[troughs])/(1/fs)/1000)                        #outwards streamwise collision velocity
Win = ((data$y[troughs]-data$y[troughs-1])/(1/fs)/1000)                         #inwards vertical collision velocity
Wout = ((data$y[troughs+1]-data$y[troughs])/(1/fs)/1000)                        #outwards streamwise collision velocity

# Store raw results in a dataframe
result = data.frame(
  saltations,                                                                                              
  t,
  Lp, 
  Hp,
  Up,
  Uin = Uin[-length(Uin)],
  Uout = Uout[-length(Uout)],
  Win = abs(Win[-length(Win)]),
  Wout = Wout[-length(Wout)]
)

# bind collision angles to raw results
results<-cbind(result, t(abs(angles[,-ncol(angles)])))

#post processing################################################################
# modes of transport threshold
result$event = "salt"                                                           #mark saltation events
result$event[abs(reld) < susthres] = "sus"                                      #mark suspension events (relative drop length < threshold)
result$event[Z0 < rollthreshp] = "roll"                                         #mark rolling events (relative height of particle < threshold)
result                                                                          #results for all events (rolling saltation and suspension)

#Final results##################################################################
#result just saltation
result[result$event=="salt",]                                                   #results for all saltation events

#Plot trajectory###############################################################
plot(data, type = 'o',xlab="x (mm)", ylab="y (mm)", xlim = c(0, 500))           #plot trajectory               
par(new = T)                                                                   
abline(v = data$x[peaks], lty = 3, col = 'blue')                                # add peaks in blue
abline(v = data$x[troughs], lty = 3, col = 'red')                               # add troughs in red  

