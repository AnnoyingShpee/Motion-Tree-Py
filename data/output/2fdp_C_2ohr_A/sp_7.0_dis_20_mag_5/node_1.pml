load 2fdp_C_2ohr_A.pdb, node_1
select region0, node_1 and resi 49-157
select region0, region0 + (node_1 and resi 169-311)
select region0, region0 + (node_1 and resi 315-385)
set_color colour0 = [0  ,0  ,255]
color colour0, region0
deselect
select region1, node_1 and resi 46-48
set_color colour1 = [255,0  ,0  ]
color colour1, region1
deselect
