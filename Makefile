DATA_DIR=./data
DOCKER_IMAGE=chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
S3_BUCKET=s3://orthrus-meta-1b449aad01564a89-inputs
ENDPOINT=--endpoint=http://storage.googleapis.com

download_raw_data:
	# mkdir -p ./data/PXD038906 # modern human coprolite spectra && \
	mkdir -p ./data/PXD027613 # ancient human coprolite spectra && \
	# conda run -n orthrus pridepy download-all-public-raw-files -a PXD038906 -o ./data/PXD038906 -p aspera && \
	conda run -n orthrus pridepy download-all-public-raw-files -a PXD027613 -o ./data/PXD027613 -p aspera

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

convert_raw_data_to_mzml_single_dir:
	@subdir=/home/j-vangoey/code/orthrus_metaproteomics/data/PXD038906; \
	if [ -d $$subdir ]; then \
		mzml_dir=$$subdir/mzML; \
		mkdir -p $$mzml_dir; \
		for file in $$subdir/*.raw; do \
			if [ -f $$file ]; then \
				echo "Processing $$file into $$mzml_dir..."; \
				docker run -it --rm -e WINEDEBUG=-all -v $$subdir:/data $(DOCKER_IMAGE) wine msconvert /data/$$(basename $$file) --outdir /data/mzML; \
			fi; \
		done; \
	fi

# Pattern rule to upload each mzML file
# find data -type f -name "*.mzML" searches for all .mzML files under the data directory.
upload_mzml_to_bucket: $(shell find data -type f -name "*.mzML" | sed 's/^/upload\//')

# Target for each mzML file
# The dirname command extracts the folder structure (PXD027613/mzML), and cut keeps only the relevant part (PXD027613).
upload/%:
	aws s3 cp $(@:upload/%=%) $(S3_BUCKET)/$(shell dirname $(@:upload/%=%) | cut -d/ -f2) $(ENDPOINT)

download_mzml_to_local:
	mkdir -p ./data/PXD027613/mzML && \
	aws s3 cp s3://orthrus-meta-1b449aad01564a89-inputs/PXD027613/20210408-Paleofeces-2604-01.mzML ./data/PXD027613/mzML

run_part1:
	conda run -n orthrus python orthrus_v1/annotated_orthrus_pt1.py
