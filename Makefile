DATA_DIR=./data
DOCKER_IMAGE=chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
S3_BUCKET=s3://orthrus-meta-1b449aad01564a89-inputs
ENDPOINT=--endpoint=http://storage.googleapis.com

download_raw_data:
	mkdir -p ./data/PXD038906 # modern human coprolite spectra && \
	mkdir -p ./data/PXD027613 # ancient human coprolite spectra && \
	conda run -n orthrus --live-stream pridepy download-all-public-raw-files -a PXD038906 -o ./data/PXD038906 -p aspera && \
	conda run -n orthrus --live-stream pridepy download-all-public-raw-files -a PXD027613 -o ./data/PXD027613 -p aspera

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
	@for subdir in data/PXD027613; do \
		mzml_dir=$$subdir/mzML; \
		if [ -d $$mzml_dir ]; then \
			folder_name=$$(basename $$subdir); \
			for file in $$mzml_dir/*.mzML; do \
				if [ -f $$file ]; then \
					base_name=$$(basename $$file); \
					echo "Uploading $$file to s3://orthrus-meta-1b449aad01564a89-inputs/$$folder_name/$$base_name..."; \
					aws s3 cp $$file s3://orthrus-meta-1b449aad01564a89-inputs/$$folder_name/$$base_name --endpoint=http://storage.googleapis.com; \
				fi; \
			done; \
		else \
			echo "Directory $$mzml_dir does not exist."; \
		fi; \
	done

download_mzml_to_local:
	mkdir -p ./data/PXD027613/mzML && \
	aws s3 cp $(S3_BUCKET)/PXD027613/20210408-Paleofeces-2604-01.mzML ./data/PXD027613/mzML && \
	aws s3 cp $(S3_BUCKET)/PXD027613/20210408-Paleofeces-2610-01.mzML ./data/PXD027613/mzML && \
	aws s3 cp $(S3_BUCKET)/PXD027613/20210408-Paleofeces-2611-01.mzML ./data/PXD027613/mzML && \
	aws s3 cp $(S3_BUCKET)/PXD027613/20210408-Paleofeces-2612-01.mzML ./data/PXD027613/mzML && \
	aws s3 cp $(S3_BUCKET)/PXD027613/20210408-Paleofeces-ControlBlank-01.mzML ./data/PXD027613/mzML

run_part1:
	conda run -n orthrus --live-stream python orthrus_v1/annotated_orthrus_pt1.py
