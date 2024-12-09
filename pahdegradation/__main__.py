import argparse

from pahdegradation.run import run

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Run PAH degradation models.")
    parser.add_argument(
        "--steps", type=int, default=50, help="Number of steps for the simulation."
    )
    parser.add_argument(
        "--optimization", action="store_true", help="Enable parameter optimization."
    )
    args = parser.parse_args()

    # Call the run function with arguments from flags
    run(steps=args.steps, optimization=args.optimization)

if __name__ == "__main__":
    main()