DOCKER_IMAGE=chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
S3_INPUT_BUCKET=s3://orthrus-meta-1b449aad01564a89-inputs
S3_OUTPUT_BUCKET=s3://orthrus-meta-1b449aad01564a89-outputs/output/b9bc0726-9636-4337-be0b-8a6dd3635247
S3_ENDPOINT=http://storage.googleapis.com
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

upload_mzml_to_s3:
	@for subdir in data/$(PRIDE_ID); do \
		mzml_dir=$$subdir/mzML; \
		if [ -d $$mzml_dir ]; then \
			folder_name=$$(basename $$subdir); \
			for file in $$mzml_dir/*.mzML; do \
				if [ -f $$file ]; then \
					base_name=$$(basename $$file); \
					echo "Uploading $$file to $(S3_INPUT_BUCKET)/$$folder_name/$$base_name..."; \
					aws s3 cp $$file $(S3_INPUT_BUCKET)/$$folder_name/$$base_name --endpoint=$(S3_ENDPOINT); \
				fi; \
			done; \
		else \
			echo "Directory $$mzml_dir does not exist."; \
		fi; \
	done

download_mzml_to_local:
	mkdir -p $(MZML_DIR) && \
	aws s3 cp $(S3_INPUT_BUCKET)/$(PRIDE_ID)  $(MZML_DIR) --recursive --endpoint=$(S3_ENDPOINT)

download_instanovo_predictions_to_local:
	mkdir -p $(MZML_DIR) && \
	aws s3 cp $(S3_OUTPUT_BUCKET)/data/$(PRIDE_ID)/mzML/ $(MZML_DIR) --recursive --endpoint=$(S3_ENDPOINT)

download_matched_fasta_to_local:
	mkdir -p $(MZML_DIR)
	aws s3 cp --recursive \
		$(S3_OUTPUT_BUCKET) \
		$(MZML_DIR) \
		--endpoint=$(S3_ENDPOINT) \
		--exclude "*" \
		--include "*_matched.fasta"

run_part1:
	conda run -n orthrus --live-stream python orthrus_v1/annotated_orthrus_pt1.py

run_part2:
	conda run -n orthrus --live-stream python orthrus_v1/annotated_orthrus_pt2.py
