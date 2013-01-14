set_color dblue, [30,144,255]
set_color fgreen, [34,139,34]
set_color fire, [178,34,34]
set_color orange, [255,165,0]
set_color orchid, [153,50,204]
set_color vred, [208,32,144]
fetch 3pbl, type=pdb1, async=0
split_state 3pbl
show cartoon
hide lines
remove solvent
util.cbc
show spheres, organic
util.cbag organic
select rare,resid 17+329+150
select common,resid 17+9+9
color orchid, rare
show sticks, rare
color orange, common
show sticks, common
