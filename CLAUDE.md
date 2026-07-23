My friend wrote some code to create a star map for a space nation RP in star-map.py. I want to improve upon it.

The game Elite: Dangerous uses a very unique algorithm named Stellar Forge to create its full-sized Milky Way galaxy map. The end goal of this project is to implement a scaled-down version of Stellar Forge that generates a new and improved star map. Stellar Forge works by splitting the galaxy into boxels, each of which has a set of astrophysical parameters such as age, distance from the galactic core, and overall mass. It then runs a scientifically-accurate simulation to distribute stars in each boxel, assigning astrophysical parameters and a random seed to each star. When a player enters a star system, the random seed and astrophysical parameters are used to deterministically generate a realistic star system.

The gif file in this directory contains an image of the old galaxy map for the RP. I would like the newly generated galaxy map to be roughly the same shape as the galaxy in the old map. Each white dot on the map is a "sector" that consists of 10,000 stars. I believe the best approach is to use the sectors as the highest-level division in the galaxy in order to broadly organize where stars are meant to be, assigning stellar formation parameters to each sector, with a separate function allowing on-demand generation of individual systems in a given sector. I want the functionality of this new map generation software to at least perform all the basic features of star-map.py. 

You will likely need:

- An "Elite: Dangerous Expert" subagent which finds information on the internet about Stellar Forge and Elite: Dangerous' galaxy generation system
- A "Coding" subagent that handles the relevant software implementation
- An "Adaptation and Verification" subagent which ensures that the outputted code produces an output that resembles the old map and matches this spec
- At least one "Astrophysicist" subagent which searches publicly available astrophysics scholarship for relevant information when asked by the other agents.

You will take the role of a coordiantor or orchestrator of these agents, ensuring they work together as efficiently as possible. If there is a more effective distribution of agents, tell me about it, tell me why it's more efficient, and then implement it. Have each subagent use a separate log file to log their work and thinking. Additionally, if there is a less token-intensive way to achieve a similarly effective result, do it. 