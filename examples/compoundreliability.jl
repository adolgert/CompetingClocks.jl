# Compound Reliability
#
# Consider a computer with parts that may fail.
#
# * CPU
# * storage drive
# * power supply
# * graphics card
# * memory module
# * fan
# * motherboard
#
# Each part fails at a different rate. Let's suppose we
# have data from a population of computers and that we've figured
# out the rate of failure of each part using survival analysis.
# In lieu of that, we can use the mean time between failure (MTBF)
# from the manufacturer.
#
# If we had a server farm full of these computers, we might ask
# some practical questions.
#
# 1. Does the model, with all of its parts, match our overall observed rate of failure?
#
# 1. When do we decide to replace a machine entirely? Do we do this pre-emptively or after
#    the n-th failure or after a particular part fails?
#
# 1. Would it make sense to always replace a few parts whenever one fails?
#
# 1. How much money would we make by spending more for a more reliable part?
#
# These questions rely on a cost model.
#
# * The total time running yields $ 0.10 per hour.
# * There is a different cost for each part.
# * There is a fixed cost for each repair.
#
