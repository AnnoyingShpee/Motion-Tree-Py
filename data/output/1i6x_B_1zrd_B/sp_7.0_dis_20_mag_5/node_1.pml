load 1i6x_B_1zrd_B.pdb, node_1
select region0, node_1 and resi 8-135
set_color colour0 = [0  ,0  ,255]
color colour0, region0
deselect
select region1, node_1 and resi 136-207
set_color colour1 = [255,0  ,0  ]
color colour1, region1
deselect
