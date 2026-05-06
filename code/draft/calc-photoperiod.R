### Photoperiod

## We're using methods from BreakingLeafBud Paper (not published) 
## Methods

# Photoperiod is calced by latitude, dayOfYear of obs, and Solar Declination
delta <- 0.409*sin(2*pi(flowering_doy - 81)/365)

# daylength (hours of daylight) is derived from the sunset hour angle (w)
# phi is lat of location
w <- acos(-tan(phi)*tan(delta))
daylength <- 24*w/pi


