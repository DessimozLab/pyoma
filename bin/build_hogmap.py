import pyoma.browser.hogidmap
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="find related hogs in new version")
    parser.add_argument("--target", required=True, help="Path to the target database")
    parser.add_argument(
        "--old", nargs="+", help="Path to databases that should be mapped"
    )
    conf = parser.parse_args()
    print(conf)

    pyoma.browser.hogidmap.build_lookup(conf.target, conf.old)
    print("done... bye bye")
