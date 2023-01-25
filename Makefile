configure:
	conda env create -f conda_env.yml
	docker pull fggutierrez2018/moma2:latest
	mkdir -p results
	mkdir -p data/pdbs
