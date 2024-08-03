--CREATE SCHEMA motion_tree_py;
SET SEARCH_PATH TO motion_tree_py, public;


CREATE TABLE IF NOT EXISTS proteins (
	protein_id		CHAR(4) CHECK (protein_id <> '    '),
	protein_chain	CHAR(1) CHECK (protein_chain <> ' '),
	protein_name	VARCHAR(100),
	dist_mat		BYTEA,
	PRIMARY KEY (protein_id, protein_chain)
);

CREATE TABLE IF NOT EXISTS motion_trees (
	protein_1		    CHAR(4),
	chain_1			    CHAR(1),
	protein_2		    CHAR(4),
	chain_2			    CHAR(1),
	spatial_proximity	DECIMAL,
	dissimilarity	    SMALLINT,
	diff_dist_mat	    BYTEA,
	link_mat            BYTEA,
	PRIMARY KEY (protein_1, chain_1, protein_2, chain_2, spatial_proximity, dissimilarity)
);

CREATE TABLE IF NOT EXISTS nodes (
	protein_1		    CHAR(4),
	chain_1			    CHAR(1),
	protein_2		    CHAR(4),
	chain_2			    CHAR(1),
	spatial_proximity	DECIMAL,
	dissimilarity	    SMALLINT,
	magnitude		    SMALLINT,
	node                SMALLINT,
	distance			DECIMAL,
	large_domain        BYTEA,
	small_domain        BYTEA,
	PRIMARY KEY (protein_1, chain_1, protein_2, chain_2, spatial_proximity, dissimilarity, magnitude, node),
	FOREIGN KEY (protein_1, chain_1, protein_2, chain_2, spatial_proximity, dissimilarity)
	REFERENCES motion_tree (protein_1, chain_1, protein_2, chain_2, spatial_proximity, dissimilarity)
	ON UPDATE CASCADE ON DELETE CASCADE
);
