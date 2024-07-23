load 4ake_A_2eck_B.pdb, node_2
select region3, node_2 and resi 1-6
select region3, region3 + (node_2 and resi 13-29)
select region3, region3 + (node_2 and resi 62-111)
select region3, region3 + (node_2 and resi 174-214)
set_color colour3 = [0  ,0  ,255]
color colour3, region3
deselect
select region4, node_2 and resi 30-61
set_color colour4 = [255,0  ,0  ]
color colour4, region4
deselect
select region5, node_2 and resi 7-12
select region5, region5 + (node_2 and resi 112-173)
set_color colour5 = [128,128,128]
color colour5, region5
deselect
