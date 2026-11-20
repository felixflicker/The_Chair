# The Chair

## Enumeration 

Enumeration of possible configurations of Chair44 and Chair611, strongly aperiodic 3D monotiles, to accompany "Matching Rules for a Three-Dimensional Strongly Aperiodic Monotile", https://arxiv.org/abs/2609.23783.

Usage: 

python3 enumerate_figs.py [--all]

Takes around 5 minutes without / 15 minutes with the --all flag. With the flag, all configurations, legal and illegal, are shown. Without, only the legal configurations are shown.

Generates gallery.html (stored here under /docs).

You can access the interactive gallery directly here:

https://felixflicker.github.io/The_Chair/gallery.html

## Growth

Three models of growth of Chair611. 

Usage: 

python3 growth.py MODE [key=value]

MODE = base, skeleton, anneal.

base: seed the base of a level-L superchair, then add a Chair611 wherever it is the only legal chair covering an empty cube next to the cluster. The result regrows the entire level-L superchair. Sample output:

https://felixflicker.github.io/The_Chair/base_L2z_f1_s1.html

skeleton: Homogeneous 3D nucleation. A candidate must be legal and must leave every empty cube it touches coverable by some other legal chair. A cube with exactly one candidate is 'forced'; if no cube is forced, choose a random candidate. After a small number of initial choices, the remaining structure, an infinite periodic 'skeleton', is forced indefinitely. Toggle to see missing Chairs. Sample output:

https://felixflicker.github.io/The_Chair/skeleton_n30000_s1.html

anneal: divide each of the 24 orientations of Chair611 into 7 cubes. Colour internal cube faces so that each pair of meeting faces gets a unique colour. This defines a set of 168 'Wang cubes'. Usually these would obey deterministic rules. Instead the simulation fixes a correct layer of cubes around a level-L Superchair, randomises every cube inside this boundary, then cools using heat-bath Monte Carlo, assigning an energy gain whenever cubes match face colours with their neighbours. The perfect level-L Superchair reforms (tested for L=1,2,3). Sample output:

https://felixflicker.github.io/The_Chair/anneal_L2_s1.html

## Chairtris

A few people have 3D printed the Chair tiling based on the Chair611 rules. Some asked if it's possible to assemble the tiling one tile at a time, or whether you need to piece together larger blocks first. It is possible to piece the tiling together one Chair at a time, while only ever dropping chairs in at right angles to the surface. But you have to get the order correct. To demonstrate this, you can play 'Chairtris': a 3D tetris-like game that teaches you how to assemble The Chair tiling. Use Q and A to rotate the incoming tile. Be sure to match black to white, and blue to blue. You may need to use 'p' to pause!

https://felixflicker.github.io/The_Chair/chairtris_standalone.html

Code written by Anthropic's Claude Fable 5.1 and 5.2.
