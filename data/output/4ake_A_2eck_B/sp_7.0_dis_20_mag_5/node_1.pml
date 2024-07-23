load 4ake_A_2eck_B.pdb, node_1
select region0, node_1 and resi 1-6
select region0, region0 + (node_1 and resi 13-111)
select region0, region0 + (node_1 and resi 174-214)
set_color colour0 = [0  ,0  ,255]
color colour0, region0
deselect
select region1, node_1 and resi 7-12
select region1, region1 + (node_1 and resi 112-173)
set_color colour1 = [255,0  ,0  ]
color colour1, region1
deselect
