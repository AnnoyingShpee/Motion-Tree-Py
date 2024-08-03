load 1so7_A_1vcu_A.pdb, node_1
select region0, node_1 and resi 4-41
select region0, region0 + (node_1 and resi 49-107)
select region0, region0 + (node_1 and resi 117-225)
select region0, region0 + (node_1 and resi 229-283)
select region0, region0 + (node_1 and resi 288-377)
set_color colour0 = [0  ,0  ,255]
color colour0, region0
deselect
select region1, node_1 and resi 108-116
set_color colour1 = [255,0  ,0  ]
color colour1, region1
deselect
