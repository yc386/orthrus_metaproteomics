import os
import s3fs
import logging
import sys

os.environ["S3FS_LOGGING_LEVEL"] = "DEBUG"
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)


def upload_to_bucket(output_path):
    """Upload results to a bucket."""
    # Only applicable when running on on https://aichor.ai/
    if "AICHOR_OUTPUT_PATH" in os.environ:
        s3_endpoint = "https://storage.googleapis.com"  # os.environ["S3_ENDPOINT"]
        s3_key = os.environ["AWS_ACCESS_KEY_ID"]
        s3_secret_key = os.environ["AWS_SECRET_ACCESS_KEY"]
        s3 = s3fs.S3FileSystem(
            client_kwargs={"endpoint_url": s3_endpoint},
            key=s3_key,
            secret=s3_secret_key,
        )
        bucket_path = f"{os.environ['AICHOR_OUTPUT_PATH']}{output_path}"
        with open(output_path, "r") as local_file, s3.open(
            bucket_path, mode="w"
        ) as bucket_file:
            bucket_file.write(local_file.read())
        print(f" 🪣 Results uploaded to {bucket_path}")
    else:
        print(" 🪣 Results not uploaded. Not running on AIchor.")


def write_large_file(filename, num_lines=1_000_000):
    """Writes a file of approximately the given target size."""
    with open(filename, "w") as file:
        for i in range(1, num_lines + 1):
            data = "*" * 100
            line = f"{i},{data}\n"
            file.write(line)

    file_size_bytes = os.path.getsize(filename)
    file_size_gb = file_size_bytes / (1024**3)
    print(f"File '{filename}' written with size: {file_size_gb:.2f} GB")


for i in range(1, 6):
    num_lines = 10**i
    write_large_file(f"large_file_{num_lines}.txt", num_lines=num_lines)
    upload_to_bucket(f"large_file_{num_lines}.txt")
