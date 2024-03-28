# Rules and Guidelines

Even though this library is just handling a bag of event clocks, there are some rules for how to call the interface so that those clocks remain consistent. There are rules that keep the simulation running, and there are rules that guarantee repeatable simulations that have the correct statistical likelihood.

 * Some simulations treat events as maybe changing state, maybe not changing state.
 * Some simulations think of an event as probabilistically changing state.
 * For a simulation in continuous time, only one event can happen at any time.
