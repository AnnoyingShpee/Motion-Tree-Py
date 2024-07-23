load 1s3i_A_2bw0_A.pdb, node_1
select region0, node_1 and resi 106-132
select region0, region0 + (node_1 and resi 136-307)
set_color colour0 = [0  ,0  ,255]
color colour0, region0
deselect
select region1, node_1 and resi 1-105
select region1, region1 + (node_1 and resi 133-135)
set_color colour1 = [255,0  ,0  ]
color colour1, region1
deselect
