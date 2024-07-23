load 4ake_A_2eck_A.pdb, node_3
select region6, node_3 and resi 30-59
set_color colour6 = [0  ,0  ,255]
color colour6, region6
deselect
select region7, node_3 and resi 60-80
set_color colour7 = [255,0  ,0  ]
color colour7, region7
deselect
select region8, node_3 and resi 1-29
select region8, region8 + (node_3 and resi 81-214)
set_color colour8 = [128,128,128]
color colour8, region8
deselect
