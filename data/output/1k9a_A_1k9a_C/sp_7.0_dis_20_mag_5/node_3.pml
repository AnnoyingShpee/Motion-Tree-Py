load 1k9a_A_1k9a_C.pdb, node_3
select region6, node_3 and resi 6-79
select region6, region6 + (node_3 and resi 175-200)
select region6, region6 + (node_3 and resi 208-241)
select region6, region6 + (node_3 and resi 251-271)
select region6, region6 + (node_3 and resi 323-324)
set_color colour6 = [0  ,0  ,255]
color colour6, region6
deselect
select region7, node_3 and resi 201-207
set_color colour7 = [255,0  ,0  ]
color colour7, region7
deselect
select region8, node_3 and resi 80-174
select region8, region8 + (node_3 and resi 242-250)
select region8, region8 + (node_3 and resi 272-322)
select region8, region8 + (node_3 and resi 325-340)
select region8, region8 + (node_3 and resi 348-450)
set_color colour8 = [128,128,128]
color colour8, region8
deselect
