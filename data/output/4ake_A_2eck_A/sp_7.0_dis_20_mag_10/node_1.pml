load 4ake_A_2eck_A.pdb, node_1
select region0, node_1 and resi 1-116
select region0, region0 + (node_1 and resi 165-214)
set_color colour0 = [0  ,0  ,255]
color colour0, region0
deselect
select region1, node_1 and resi 117-164
set_color colour1 = [255,0  ,0  ]
color colour1, region1
deselect
