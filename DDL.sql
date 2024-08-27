--CREATE SCHEMA motion_tree_py;
SET SEARCH_PATH TO motion_tree_py, public;


CREATE TABLE IF NOT EXISTS proteins (
	protein_1		VARCHAR(15),
	chain_1			CHAR(1),
	protein_2		VARCHAR(15),
	chain_2			CHAR(1),
	protein_name    VARCHAR(100),
	rmsd            DECIMAL,
	diff_dist_mat   BYTEA,
	PRIMARY KEY (protein_1, chain_1, protein_2, chain_2)
);

CREATE TABLE IF NOT EXISTS motion_trees (
	protein_1		    VARCHAR(15),
	chain_1			    CHAR(1),
	protein_2		    VARCHAR(15),
	chain_2			    CHAR(1),
	spatial_proximity	DECIMAL,
	sequence_identity   DECIMAL,
	time_taken          DECIMAL,
	link_mat            BYTEA,
	PRIMARY KEY (protein_1, chain_1, protein_2, chain_2, spatial_proximity),
	FOREIGN KEY (protein_1, chain_1, protein_2, chain_2)
	REFERENCES proteins (protein_1, chain_1, protein_2, chain_2)
	ON UPDATE CASCADE ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS nodes (
	protein_1		    VARCHAR(15),
	chain_1			    CHAR(1),
	protein_2		    VARCHAR(15),
	chain_2			    CHAR(1),
	spatial_proximity	DECIMAL,
	small_node_size     SMALLINT,
	cluster_size	    SMALLINT,
	magnitude		    SMALLINT,
	node                SMALLINT,
	distance			DECIMAL,
	large_domain        BYTEA,
	small_domain        BYTEA,
	PRIMARY KEY (protein_1, chain_1, protein_2, chain_2, spatial_proximity, small_node_size, cluster_size, magnitude, node),
	FOREIGN KEY (protein_1, chain_1, protein_2, chain_2, spatial_proximity)
	REFERENCES motion_tree (protein_1, chain_1, protein_2, chain_2, spatial_proximity)
	ON UPDATE CASCADE ON DELETE CASCADE
);
