# # Reliability of a Group
#
# Generalized semi-Markov processes are often used for reliability modeling.
# A simple reliability model for some vehicle might say it is either
# working or in repair. There is a distribution of the working time and
# a distribution of the repair time. Some models make repairs take A
# fixed amount of time. The Fleck library can certainly model an individual,
# but let's consider a group of vehicles.
#
# We need to use ten vehicles a day, and we start with eighteen vehicles.
# As vehicles go into repair, the ten for the day are chosen from the
# remaining usable vehicles.
#

# The simplest model for one vehicle:
#
# * Common choices for the failure distribution are Weibull and Lognormal.
# * Repair distribution is a fixed time, so it's dirac-distributed.
#
# How do we account for how much the vehicle is used?
# We could make some kind of age function which depends on the
# number of days it is used and the amount of time each day.
# In the language of the GSMP, every day the vehicle is used,
# we will enable the failure transition. When we disable it,
# we will ensure that we remember the total length of time that
# failure transition has run so far before failing.
#
