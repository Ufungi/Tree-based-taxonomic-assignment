import json
import argparse

def clean_jplace(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    data["placements"] = [p for p in data["placements"] if not any("nan" in str(v) for lst in p["p"] for v in lst)]

    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove placements containing NaN from query.jplace")
    parser.add_argument("input", help="Path to the input query.jplace file")
    parser.add_argument("output", help="Path to save the cleaned query.jplace file")

    args = parser.parse_args()
    clean_jplace(args.input, args.output)
