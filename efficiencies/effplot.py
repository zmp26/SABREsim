
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_columns.py <path_to_txt_file>")
        sys.exit(1)

    file_path = sys.argv[1]

    try:
        # Load data, skip header row
        data = np.loadtxt(file_path, skiprows=1)

        # Column indices are 0-based
        x = data[:, 0]   # 1st column
        y = data[:, 7]   # 8th column

        plt.figure(figsize=(8,4))
        plt.plot(x, y)
        plt.xlabel(r"$^{9}$B Excitation Energy (keV)")
        plt.ylabel("3-Particle Detection Efficiency")
        plt.grid(True)
        plt.show()

    except Exception as e:
        print(f"Error reading or plotting file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()