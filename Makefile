DATA_DIR=./data
DOCKER_IMAGE=chambm/pwiz-skyline-i-agree-to-the-vendor-licenses

download_data:
	# mkdir -p ./data/PXD038906 # modern human coprolite spectra && \
	mkdir -p ./data/PXD027613 # ancient human coprolite spectra && \
	# conda run -n orthrus pridepy download-all-public-raw-files -a PXD038906 -o ./data/PXD038906 -p aspera && \
	conda run -n orthrus pridepy download-all-public-raw-files -a PXD027613 -o ./data/PXD027613 -p aspera

convert_data:
	@for subdir in $(DATA_DIR)/*; do \
		if [ -d $$subdir ]; then \
			mzml_dir=$$subdir/mzML; \
			mkdir -p $$mzml_dir; \
			for file in $$subdir/*.raw; do \
				if [ -f $$file ]; then \
					echo "Processing $$file into $$mzml_dir..."; \
					docker run -it --rm -e WINEDEBUG=-all -v $$subdir:/data $(DOCKER_IMAGE) wine msconvert /data/$$(basename $$file) --outdir /data/mzML; \
				fi; \
			done; \
		fi; \
	done

run_part1:
	conda run -n orthrus python orthrus_v1/annotated_orthrus_pt1.py
