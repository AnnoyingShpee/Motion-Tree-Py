load 1k9a_A_1k9a_C.pdb, node_4
select region9, node_4 and resi 6-79
select region9, region9 + (node_4 and resi 175-186)
select region9, region9 + (node_4 and resi 229-241)
set_color colour9 = [0  ,0  ,255]
color colour9, region9
deselect
select region10, node_4 and resi 187-200
select region10, region10 + (node_4 and resi 208-228)
select region10, region10 + (node_4 and resi 251-271)
select region10, region10 + (node_4 and resi 323-324)
set_color colour10 = [255,0  ,0  ]
color colour10, region10
deselect
select region11, node_4 and resi 80-174
select region11, region11 + (node_4 and resi 201-207)
select region11, region11 + (node_4 and resi 242-250)
select region11, region11 + (node_4 and resi 272-322)
select region11, region11 + (node_4 and resi 325-340)
select region11, region11 + (node_4 and resi 348-450)
set_color colour11 = [128,128,128]
color colour11, region11
deselect
