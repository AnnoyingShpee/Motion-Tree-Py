load 4qfk_A_4qfl_A.pdb, node_1
select region0, node_1 and resi 30-289
select region0, region0 + (node_1 and resi 507-534)
set_color colour0 = [0  ,0  ,255]
color colour0, region0
deselect
select region1, node_1 and resi 290-506
set_color colour1 = [255,0  ,0  ]
color colour1, region1
deselect
