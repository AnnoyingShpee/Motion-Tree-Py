load 1k9a_A_1k9a_C.pdb, node_2
select region3, node_2 and resi 242-250
select region3, region3 + (node_2 and resi 272-322)
select region3, region3 + (node_2 and resi 325-340)
select region3, region3 + (node_2 and resi 348-450)
set_color colour3 = [0  ,0  ,255]
color colour3, region3
deselectselect region4, node_2 and resi 6-79
select region4, region4 + (node_2 and resi 175-241)
select region4, region4 + (node_2 and resi 251-271)
select region4, region4 + (node_2 and resi 323-324)
set_color colour4 = [255,0  ,0  ]
color colour4, region4
deselectselect region5, node_2 and resi 80-174
set_color colour5 = [128,128,128]
color colour5, region5
deselect