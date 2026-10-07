# The Chair

Enumeration of possible configurations of Chair44 and Chair611, strongly aperiodic 3D monotiles, to accompany "Matching Rules for a Three-Dimensional Strongly Aperiodic Monotile", https://arxiv.org/abs/2609.23783.

Usage: 

python3 enumerate_figs.py [--all]

Takes around 5 minutes without / 15 minutes with the --all flag. With the flag, all configurations, legal and illegal, are shown. Without, only the legal configurations are shown.

Generates gallery.html (stored here under /docs).

You can access the interactive gallery directly here:

https://felixflicker.github.io/The_Chair/gallery.html

A few people have 3D printed the Chair tiling based on the Chair611 rules. Some asked if it's possible to assemble the tiling one tile at a time, or whether you need to piece together larger blocks first. It is possible to piece the tiling together one Chair at a time, while only ever dropping chairs in at right angles to the surface. But you have to get the order correct. To demonstrate this, you can play 'Chairtris': a 3D tetris-like game that teaches you how to assemble The Chair tiling. Use Q and A to rotate the incoming tile. Be sure to match black to white, and blue to blue. You may need to use 'p' to pause!

https://felixflicker.github.io/The_Chair/chairtris_standalone.html

Code written by Anthropic's Claude Fable 5.1.
