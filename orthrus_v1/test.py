import os
import s3fs


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


with open("output.txt", "w") as f:
    f.write("Hello, world!")

upload_to_bucket("output.txt")
