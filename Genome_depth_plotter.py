## Author: Wentao Chen
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import subprocess
from matplotlib.colors import to_rgb, LinearSegmentedColormap

def calculate_coverage_depth(bam_file):
    try:
        print(f"Processing BAM file: {bam_file}")

        # Execute samtools depth to get coverage information
        cmd = ["samtools", "depth", bam_file]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        out, err = process.communicate()

        if process.returncode != 0:
            print(f"Error in samtools depth: {err}")
            return np.array([]), 0

        # Parse the output and create coverage array
        coverage_data = [line.split('\t') for line in out.splitlines()]
        if not coverage_data:
            print(f"No coverage data found in {bam_file}")
            return np.array([]), 0

        # Find the max position as reference length
        ref_length = max(int(line.split('\t')[1]) for line in out.splitlines())
        coverage = np.zeros(ref_length, dtype=np.int32)
        for line in coverage_data:
            pos = int(line[1])
            depth = int(line[2])
            coverage[pos - 1] = depth

        print(f"Coverage calculation completed for {bam_file}.")
        return coverage, ref_length

    except Exception as e:
        print(f"Error processing BAM file {bam_file}: {e}")
        return np.array([]), 0

def generate_color_gradient(center_color, num_colors):
    """
    Generate a color gradient based on a center color.
    :param center_color: Central color for the gradient
    :param num_colors: Number of colors to generate
    :return: Array of colors
    """
    print("Generating color gradient...")
    base_color = to_rgb(center_color)
    cmap = LinearSegmentedColormap.from_list("custom_cmap", [(1, 1, 1), base_color, (0, 0, 0)])
    return [cmap(i) for i in np.linspace(0.2, 0.8, num_colors)]

def set_x_axis(ax, ref_length):
    """
    Set up the X-axis of the plot.
    :param ax: Axis object for the plot
    :param ref_length: Length of the reference genome
    """
    ax.set_xlim(1, ref_length)
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useOffset=False))
    ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

def plot_genome_depth(coverage, ref_length, adjusted_font_size_x, adjusted_font_size_y, color, alpha, max_depth, depth_indicator, ax, show_title, grid, show_x_label):
    """
    Plot the genome depth.
    :param coverage: Coverage array
    :param ref_length: Reference genome length
    :param adjusted_font_size_x: Font size for the X-axis label
    :param adjusted_font_size_y: Font size for the Y-axis label
    :param color: Color for the plot
    :param alpha: Alpha (transparency) for the plot area
    :param max_depth: Maximum depth for the Y-axis
    :param depth_indicator: Depth value to indicate with a horizontal line
    :param ax: Axis object for the plot
    :param show_title: Whether to display the title
    :param grid: Whether to display the grid
    :param show_x_label: Whether to show X-axis label
    """
    print("Creating genome depth plot...")
    x = np.arange(1, len(coverage) + 1)
    set_x_axis(ax, ref_length)
    ax.set_ylim(0, max_depth)
    if show_x_label:
        ax.set_xlabel("Reference Genome Position", fontsize=adjusted_font_size_x)
    ax.set_ylabel("Mapped Read Depth (X)", fontsize=adjusted_font_size_y)
    ax.plot(x, coverage, color=color, linestyle='-', alpha=alpha,linewidth=0.1)
    ax.fill_between(x, 0, coverage, color=color, alpha=alpha)

    if depth_indicator is not None:
        ax.axhline(depth_indicator, color='black', linestyle='--',linewidth=0.5)

    if show_title:
        ax.set_title("Coverage Plot", fontsize=adjusted_font_size_x)

    if grid:
        ax.grid(color='gray', linestyle='--', linewidth=0.5)

    ax.set_facecolor('white')
    print("Genome depth plot created.")

def main():
    parser = argparse.ArgumentParser(description="Plot genome depth from BAM files.")
    parser.add_argument("--bam_files", nargs="+", required=True, help="List of paths to the BAM files.")
    parser.add_argument("--center_color", default="#0000FF", help="Central color for the plot area.")
    parser.add_argument("--area_alpha", type=float, default=0.5, help="Alpha value for the plot area.")
    parser.add_argument("--max_depth", type=int, default=100, help="Maximum depth for the y-axis.")
    parser.add_argument("--depth_indicator", type=int, help="Depth value to indicate with a horizontal line.")
    parser.add_argument("--adjusted_font_size_x", type=int, default=12, help="Font size for the X-axis label.")
    parser.add_argument("--adjusted_font_size_y", type=int, default=9, help="Font size for the Y-axis label.")
    parser.add_argument("--fig_width", type=float, default=15, help="Width of the plot.")
    parser.add_argument("--fig_height", type=float, default=3, help="Height of the plot.")
    parser.add_argument("--show_title", action='store_true', help="Whether to display the title or not.")
    parser.add_argument("--grid", action='store_true', help="Whether to display grid lines or not.")
    parser.add_argument("--plot_mode", choices=["single", "side_by_side"], default="single", help="Plot mode: single or side by side.")
    parser.add_argument("--output_file", required=True, help="Base name for output files.")
    parser.add_argument("--output_formats", nargs="+", choices=['png', 'pdf', 'tiff'], default=['tiff'], help="List of output file formats.")
    args = parser.parse_args()

    color_gradient = generate_color_gradient(args.center_color, len(args.bam_files))
    plot_size = (args.fig_width, args.fig_height)

    # Process each BAM file and create plots
    if args.plot_mode == 'single':
        for bam_file, color in zip(args.bam_files, color_gradient):
            coverage, ref_length = calculate_coverage_depth(bam_file)
            fig, ax = plt.subplots(figsize=plot_size)
            plot_genome_depth(coverage, ref_length, args.adjusted_font_size_x, args.adjusted_font_size_y, color, args.area_alpha, args.max_depth, args.depth_indicator, ax, args.show_title, args.grid, show_x_label=True)
            for ext in args.output_formats:
                output_filename = f"{args.output_file}_{os.path.splitext(os.path.basename(bam_file))[0]}.{ext}"
                plt.savefig(output_filename, bbox_inches='tight', dpi=600)
            plt.close()
            print(f"Plot generated for {bam_file} in formats: {', '.join(args.output_formats)}")

    elif args.plot_mode == 'side_by_side':
        # Creating a subplot for each BAM file
        fig, axes = plt.subplots(len(args.bam_files), 1, figsize=(plot_size[0], plot_size[1] * len(args.bam_files)), sharex=True)

        # Reducing the vertical spacing between subplots
        fig.subplots_adjust(hspace=0.1)  # Adjust the vertical space between subplots

        for i, (bam_file, ax, color) in enumerate(zip(args.bam_files, axes, color_gradient)):
            coverage, ref_length = calculate_coverage_depth(bam_file)
            show_x_label = (i == len(args.bam_files) - 1)  # Display the x-axis label only for the last subplot
            plot_genome_depth(coverage, ref_length, args.adjusted_font_size_x, args.adjusted_font_size_y, color, args.area_alpha, args.max_depth, args.depth_indicator, ax, args.show_title, args.grid, show_x_label)

        # Saving the entire figure
        for ext in args.output_formats:
            output_filename = f"{args.output_file}_side_by_side.{ext}"
            plt.savefig(output_filename, bbox_inches='tight', dpi=600)
        plt.close()
        print(f"Side by side plots generated for {args.output_file} in formats: {', '.join(args.output_formats)}")

if __name__ == "__main__":
    main()
