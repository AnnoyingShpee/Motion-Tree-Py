load 2fdp_C_2ohr_A.pdb, node_2
select region3, node_2 and resi 49-63
select region3, region3 + (node_2 and resi 76-157)
select region3, region3 + (node_2 and resi 169-311)
select region3, region3 + (node_2 and resi 315-385)
set_color colour3 = [0  ,0  ,255]
color colour3, region3
deselect
select region4, node_2 and resi 64-75
set_color colour4 = [255,0  ,0  ]
color colour4, region4
deselect
select region5, node_2 and resi 46-48
set_color colour5 = [128,128,128]
color colour5, region5
deselect
