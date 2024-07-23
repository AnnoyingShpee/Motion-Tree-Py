load 1oj7_A_1oj7_C.pdb, node_1
select region0, node_1 and resi -2-6
select region0, region0 + (node_1 and resi 12-19)
select region0, region0 + (node_1 and resi 174-387)
set_color colour0 = [0  ,0  ,255]
color colour0, region0
deselectselect region1, node_1 and resi 7-11
select region1, region1 + (node_1 and resi 20-173)
set_color colour1 = [255,0  ,0  ]
color colour1, region1
deselect