load 1so7_A_1vcu_A.pdb, node_3
select region6, node_3 and resi 4-14
select region6, region6 + (node_3 and resi 19-41)
select region6, region6 + (node_3 and resi 51-107)
select region6, region6 + (node_3 and resi 117-225)
select region6, region6 + (node_3 and resi 229-283)
select region6, region6 + (node_3 and resi 288-377)
set_color colour6 = [0  ,0  ,255]
color colour6, region6
deselect
select region7, node_3 and resi 15-18
set_color colour7 = [255,0  ,0  ]
color colour7, region7
deselect
select region8, node_3 and resi 49-50
select region8, region8 + (node_3 and resi 108-116)
set_color colour8 = [128,128,128]
color colour8, region8
deselect
