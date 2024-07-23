load 4ake_A_2eck_B.pdb, node_3
select region6, node_3 and resi 121-158
set_color colour6 = [0  ,0  ,255]
color colour6, region6
deselect
select region7, node_3 and resi 7-12
select region7, region7 + (node_3 and resi 112-120)
select region7, region7 + (node_3 and resi 159-173)
set_color colour7 = [255,0  ,0  ]
color colour7, region7
deselect
select region8, node_3 and resi 1-6
select region8, region8 + (node_3 and resi 13-111)
select region8, region8 + (node_3 and resi 174-214)
set_color colour8 = [128,128,128]
color colour8, region8
deselect
