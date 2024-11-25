DOCKER_IMAGE=chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
DATA_DIR=./data
PRIDE_ID=PXD027613# ancient human coprolite spectra
# PRIDE_ID=PXD038906 # modern human coprolite spectra
MZML_DIR=./data/$(PRIDE_ID)/mzML


download_raw_data:
	mkdir -p ./data/$(PRIDE_ID)  && \
	conda run -n orthrus --live-stream pridepy download-all-public-raw-files -a $(PRIDE_ID) -o ./data/$(PRIDE_ID) -p aspera

convert_raw_data_to_mzml:
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
	conda run -n orthrus --live-stream python orthrus_v1/annotated_orthrus_pt1.py

run_part2:
	conda run -n orthrus --live-stream python orthrus_v1/annotated_orthrus_pt2.py
