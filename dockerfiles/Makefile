docker: .htslib .augustus .busco .python3

.htslib: htslib.Dockerfile
	docker build --file $< -t darcyabjones/htslib:1.9 .
	touch $@

.augustus: augustus.Dockerfile
	docker build --file $< -t darcyabjones/augustus:latest .
	touch $@

.busco: busco.Dockerfile .augustus .python3
	docker build --file $< -t darcyabjones/busco:latest .
	touch $@

.python3: python3.Dockerfile
	docker build --file $< -t darcyabjones/python3:latest .
	touch $@


