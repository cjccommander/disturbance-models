## Calculating Elasticities ##

# Read in data [Dissertation Ch. 1]
Data = read.csv("TypeIICyclicMortSD0_thr5pct_2500sims.csv", header= FALSE)
# header=FALSE because first row does not contain variable names
# row -> magnitude, column -> frequency

Data_prey = Data[ ,c(1:11)] # prey results
#Data_pred = Data[ ,c(12:22)] # predator results

Mag = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
#Mag = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) # if dist. applied to abundance
Freq = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

#Formulas for calculating elasticities for mag. and freq:
# from row i to row i-1, and column j to column j-1 --> "DOWNHILL"
# Freq: ((E(i,j-1)-E(i,j))/E(i,j))/((Freq(j-1)-Freq(j))/Freq(j))
# Mag:  ((E(i-1,j)-E(i,j))/E(i,j))/((Mag(i-1)-Mag(i))/Mag(i))
# from row i-1 to row i, and column j-1 to column j --> "UPHILL"
# Freq: ((E(i,j)-E(i,j-1))/E(i,j-1))/((Freq(j)-Freq(j-1))/Freq(j-1))
# Mag:  ((E(i,j)-E(i-1,j))/E(i-1,j))/((Mag(i)-Mag(i-1))/Mag(i-1))

# pre-allocate
eF = matrix(NA, nrow = 11, ncol = 11) #nrow=10 if dist. applied to abundance

# for loop to calculate elasticities
for (i in 1:11) # i in 1:10 if dist. applied to abundance
  for (j in 2:11)
  {
    eF[i,j] = ((Data_prey[i,j-1]-Data_prey[i,j])/Data_prey[i,j]) / 
      ((Freq[j-1]-Freq[j])/Freq[j]) # DOWNHILL
    #eF[i,j] = ((Data_prey[i,j]-Data_prey[i,j-1])/Data_prey[i,j-1]) / 
      #((Freq[j]-Freq[j-1])/Freq[j-1]) # UPHILL
  }

eF[is.nan(eF)] = 0 #replaces NaN's (resulting from 0/0) with 0's
eF[eF < 0] = 0 # replace negative values with 0
eF = abs(eF) # this is here because I kept getting negative thetas

# make .csv file for eF (match prefix (before _eF) to file above)
#write.csv(eF, file = "_eF.csv")

# pre-allocate
eM = matrix(NA, nrow = 11, ncol = 11) #nrow=10 if dist. applied to abundance

# for loop to calculate elasticities
for (i in 2:11) # i in 2:10 if dist. applied to abundance
  for (j in 1:11)
  {
    eM[i,j] = ((Data_prey[i-1,j]-Data_prey[i,j])/Data_prey[i,j]) / 
      ((Mag[i-1]-Mag[i])/Mag[i]) # DOWNHILL
    #eM[i,j] = ((Data_prey[i,j]-Data_prey[i-1,j])/Data_prey[i-1,j]) / 
      #((Mag[i]-Mag[i-1])/Mag[i-1]) # UPHILL
  }

eM[is.nan(eM)] = 0 #replaces NaN's (resulting from Div/0) with 0's
eM[eM < 0] = 0 # replace negative values with 0
eM = abs(eM) # this is here because I kept getting negative thetas

# make .csv file for eM (match prefix (before _eM) to file above)
#write.csv(eM, file = "_eM.csv")

#plot
# make .png file for plot (match prefix (before .png) to file above)
#png(filename=".png")
plot(eF,eM,xlim=c(0,8),ylim=c(0,8))
abline(0,1) # slope of 1 indicates theta = 45 degrees
#dev.off()

# Convert Cartesian coordinates to polar:
# eF is x and eM is y; convert to get (r,theta)
# Use Pythagorean theorem to get r 
# and solve for theta using the inverse tangent function -> theta = tan^-1(y/x)
# calculate Rs and Thetas
# pre-allocate
Rs = matrix(NA, nrow = 11, ncol = 11) #nrow=10 if dist. applied to abundance
Thetas_rad = matrix(NA, nrow = 11, ncol = 11) #nrow=10 if dist. applied to abundance
Thetas = matrix(NA, nrow = 11, ncol = 11) #nrow=10 if dist. applied to abundance
  
# Use Radians to Degrees function (atan function gives radians)
rad2deg = function(rad) {
    return((180 * rad) / pi)
  }

# for loop
for (i in 1:11) # i in 1:10 if dist. applied to abundance
  for (j in 1:11)
  {
Rs[i,j] = sqrt(eF[i,j]^2 + eM[i,j]^2)
Thetas_rad[i,j] = atan(eM[i,j]/eF[i,j])
Thetas[i,j] = rad2deg(Thetas_rad[i,j])
  }

Thetas[is.nan(Thetas)] = 45 #atan(0/0) and atan(Inf/Inf) result in NaN; replaces NaN's with 45 degrees

# Note: atan(0/integer) results in theta of 0 degrees; atan(integer/0) results in theta of 90 degrees
# when y=0, vector points in x direction (0 degrees); when x=0, vector points in y direction (90 degrees)

# The angle (relative to horizontal) theta of that combined vector
# will tell you the relative elasticity (is it > 45º or < 45º?)

# make .csv file for Rs (match prefix (before _Rs) to file above)
#write.csv(Rs, file = "_Rs.csv")
# make .csv file for Thetas (match prefix (before _Thetas) to file above)
write.csv(Thetas, file = "TypeIICyclicMortSD0_thr5pct_2500sims_Thetas.csv")
