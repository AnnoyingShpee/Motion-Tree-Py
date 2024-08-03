load 1so7_A_1vcu_A.pdb, node_2
select region3, node_2 and resi 4-41
select region3, region3 + (node_2 and resi 51-107)
select region3, region3 + (node_2 and resi 117-225)
select region3, region3 + (node_2 and resi 229-283)
select region3, region3 + (node_2 and resi 288-377)
set_color colour3 = [0  ,0  ,255]
color colour3, region3
deselect
select region4, node_2 and resi 49-50
set_color colour4 = [255,0  ,0  ]
color colour4, region4
deselect
select region5, node_2 and resi 108-116
set_color colour5 = [128,128,128]
color colour5, region5
deselect
