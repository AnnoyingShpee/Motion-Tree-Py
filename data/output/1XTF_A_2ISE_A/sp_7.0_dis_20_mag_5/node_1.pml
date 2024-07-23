load 1XTF_A_2ISE_A.pdb, node_1
select region0, node_1 and resi 2-420
set_color colour0 = [0  ,0  ,255]
color colour0, region0
deselect
select region1, node_1 and resi 421-422
set_color colour1 = [255,0  ,0  ]
color colour1, region1
deselect
