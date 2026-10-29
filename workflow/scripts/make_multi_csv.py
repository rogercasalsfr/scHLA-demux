#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Write a Cell Ranger multi CSV for one GEX library and an optional TCR library."
    )
    parser.add_argument("--gex-reference", required=True)
    parser.add_argument("--gex-fastqs", required=True)
    parser.add_argument("--gex-fastq-id", required=True)
    parser.add_argument("--vdj-reference")
    parser.add_argument("--tcr-fastqs")
    parser.add_argument("--tcr-fastq-id")
    parser.add_argument("--expect-cells")
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    lines = [
        "[gene-expression]",
        f"reference,{args.gex_reference}",
        "create-bam,true",
    ]
    if args.expect_cells:
        lines.append(f"expect-cells,{args.expect_cells}")

    if args.tcr_fastqs:
        lines += ["", "[vdj]", f"reference,{args.vdj_reference}"]

    lines += [
        "",
        "[libraries]",
        "fastq_id,fastqs,feature_types",
        f"{args.gex_fastq_id},{args.gex_fastqs},gene expression",
    ]
    if args.tcr_fastqs:
        lines.append(f"{args.tcr_fastq_id},{args.tcr_fastqs},vdj-t")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")
    print(f"Wrote Cell Ranger multi CSV to {output}")


if __name__ == "__main__":
    main()
