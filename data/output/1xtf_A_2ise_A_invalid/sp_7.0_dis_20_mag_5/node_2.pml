load 1XTF_A_2ISE_A.pdb, node_2
select region3, node_2 and resi 2-245
select region3, region3 + (node_2 and resi 256-420)
set_color colour3 = [0  ,0  ,255]
color colour3, region3
deselect
select region4, node_2 and resi 246-255
set_color colour4 = [255,0  ,0  ]
color colour4, region4
deselect
select region5, node_2 and resi 421-422
set_color colour5 = [128,128,128]
color colour5, region5
deselect
