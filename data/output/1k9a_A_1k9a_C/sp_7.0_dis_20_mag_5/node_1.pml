load 1k9a_A_1k9a_C.pdb, node_1
select region0, node_1 and resi 6-79
select region0, region0 + (node_1 and resi 175-340)
select region0, region0 + (node_1 and resi 348-450)
set_color colour0 = [0  ,0  ,255]
color colour0, region0
deselectselect region1, node_1 and resi 80-174
set_color colour1 = [255,0  ,0  ]
color colour1, region1
deselect