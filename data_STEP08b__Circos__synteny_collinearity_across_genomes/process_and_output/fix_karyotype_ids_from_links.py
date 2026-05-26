#!/usr/bin/env python3

from pathlib import Path
import sys

def get_link_ids(links_path: Path):
    with links_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 4:
                return fields[0], fields[3]
    return None, None


def fix_karyotype(karyo_in: Path, karyo_out: Path, new_id: str):
    with karyo_in.open() as fin, karyo_out.open("w") as fout:
        for line in fin:
            if line.startswith("chr"):
                parts = line.rstrip("\n").split()
                # Circos karyotype format:
                # chr - <id> <label> <start> <end> <color>
                parts[2] = new_id
                fout.write(" ".join(parts) + "\n")
            else:
                fout.write(line)


def main():
    links = Path("links.txt")
    if not links.exists():
        sys.exit("ERROR: links.txt not found")

    id1, id2 = get_link_ids(links)
    if id1 is None:
        print("links.txt is empty, no changes made")
        return

    print(f"Using IDs from links.txt: {id1}, {id2}")

    fix_karyotype(
        Path("genome1.karyotype.txt"),
        Path("genome1.karyotype.fixed.txt"),
        id1
    )

    fix_karyotype(
        Path("genome2.karyotype.txt"),
        Path("genome2.karyotype.fixed.txt"),
        id2
    )

    print("Wrote:")
    print("  genome1.karyotype.fixed.txt")
    print("  genome2.karyotype.fixed.txt")


if __name__ == "__main__":
    main()