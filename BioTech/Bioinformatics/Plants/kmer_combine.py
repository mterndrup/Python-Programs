import os
import pandas as pd
import re
from collections import defaultdict
from datetime import datetime
import sys
import json


class Logger:
    def __init__(self, log_file_path):
        self.log_file = open(log_file_path, 'w', encoding='utf-8')

    def write(self, message):
        print(message, end='')
        self.log_file.write(message)
        self.log_file.flush()

    def close(self):
        self.log_file.close()


def load_cluster_colors(folder_path):
    """Load cluster colors from k6 JSON file"""
    color_file_path = os.path.join(folder_path, "cluster_colors_k6.json")
    if os.path.exists(color_file_path):
        try:
            with open(color_file_path, 'r') as f:
                color_data = json.load(f)
            return color_data.get('colors', {})
        except Exception as e:
            print(f"Warning: Could not load cluster colors from {color_file_path}: {e}")
            return {}
    else:
        print(f"Warning: Cluster color file not found: {color_file_path}")
        return {}


def extract_primer_from_filename(source_file):
    """Extract primer name from source filename"""
    # Look for common primer patterns at the end before .fastq
    if '_18S.fastq' in source_file:
        return '18S'
    elif '_FITS.fastq' in source_file:
        return 'FITS'
    elif '_PITS.fastq' in source_file:
        return 'PITS'
    elif '_ITS.fastq' in source_file:
        return 'ITS'
    elif '_COI.fastq' in source_file:
        return 'COI'
    else:
        # Try to extract the last part before .fastq
        parts = source_file.replace('.fastq', '').split('_')
        if len(parts) > 1:
            potential_primer = parts[-1]
            # Return it if it looks like a primer (all caps, 2-6 characters)
            if potential_primer.isupper() and 2 <= len(potential_primer) <= 6:
                return potential_primer
        return 'Unknown'


def analyze_primer_distribution(df):
    """Analyze primer distribution per cluster"""
    # Extract primer information
    df['primer'] = df['source_file'].apply(extract_primer_from_filename)

    # Calculate primer percentages per cluster
    cluster_primer_stats = {}

    for cluster_id in df['cluster'].unique():
        cluster_data = df[df['cluster'] == cluster_id]
        primer_counts = cluster_data['primer'].value_counts()
        total_sequences = len(cluster_data)

        primer_percentages = {}
        for primer, count in primer_counts.items():
            percentage = (count / total_sequences) * 100
            primer_percentages[primer] = {'count': count, 'percentage': percentage}

        cluster_primer_stats[cluster_id] = {
            'total_sequences': total_sequences,
            'primers': primer_percentages
        }

    return cluster_primer_stats


def rgba_to_description(rgba_values):
    """Convert RGBA values to a color description"""
    if len(rgba_values) < 3:
        return "Unknown Color"

    r, g, b = rgba_values[0], rgba_values[1], rgba_values[2]

    # Convert to 0-255 range
    r_255 = int(r * 255)
    g_255 = int(g * 255)
    b_255 = int(b * 255)

    # Simple color categorization
    if r > 0.8 and g < 0.3 and b < 0.3:
        return f"Red (RGB: {r_255}, {g_255}, {b_255})"
    elif g > 0.8 and r < 0.3 and b < 0.3:
        return f"Green (RGB: {r_255}, {g_255}, {b_255})"
    elif b > 0.8 and r < 0.3 and g < 0.3:
        return f"Blue (RGB: {r_255}, {g_255}, {b_255})"
    elif r > 0.8 and g > 0.8 and b < 0.3:
        return f"Yellow (RGB: {r_255}, {g_255}, {b_255})"
    elif r > 0.8 and b > 0.8 and g < 0.3:
        return f"Magenta (RGB: {r_255}, {g_255}, {b_255})"
    elif g > 0.8 and b > 0.8 and r < 0.3:
        return f"Cyan (RGB: {r_255}, {g_255}, {b_255})"
    elif r > 0.6 and g > 0.3 and b < 0.3:
        return f"Orange (RGB: {r_255}, {g_255}, {b_255})"
    elif r > 0.5 and g < 0.3 and b > 0.5:
        return f"Purple (RGB: {r_255}, {g_255}, {b_255})"
    elif r < 0.2 and g < 0.2 and b < 0.2:
        return f"Black (RGB: {r_255}, {g_255}, {b_255})"
    elif r > 0.8 and g > 0.8 and b > 0.8:
        return f"White (RGB: {r_255}, {g_255}, {b_255})"
    elif r > 0.4 and g > 0.4 and b > 0.4:
        return f"Gray (RGB: {r_255}, {g_255}, {b_255})"
    else:
        return f"Mixed Color (RGB: {r_255}, {g_255}, {b_255})"


def get_kmer_from_filename(filename):
    """Extract k-mer value from filename like 'barcode74_k6_clusters.csv'"""
    match = re.search(r'_k(\d+)_', filename)
    return int(match.group(1)) if match else None


def analyze_barcode_clusters():
    # Program info
    script_name = os.path.basename(__file__)
    start_time = datetime.now()

    # Create log file with timestamp
    timestamp = start_time.strftime("%Y%m%d_%H%M%S")
    log_dir = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Logs"
    os.makedirs(log_dir, exist_ok=True)
    log_file_path = os.path.join(log_dir, f"barcode_cluster_analysis_{timestamp}.log")

    # Initialize logger
    logger = Logger(log_file_path)

    # Program header
    header = f"""
{'=' * 60}
Script: {script_name}
{'=' * 60}
Start Time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}
Log File: {log_file_path}
{'=' * 60}

"""
    logger.write(header)

    try:
        # Get barcode number from user
        barcode_input = input("Enter barcode number (e.g., 01-96): ")

        # Construct folder path
        base_path = r"C:\Users\ketgl\OneDrive\Desktop\Sandbox\Cactus"
        folder_path = os.path.join(base_path, f"barcode{barcode_input}")

        if not os.path.exists(folder_path):
            error_msg = f"Error: Folder {folder_path} does not exist!\n"
            logger.write(error_msg)
            return

        logger.write(f"Analyzing clusters in: {folder_path}\n")
        logger.write("-" * 50 + "\n")

        # Load cluster colors from k6 JSON file
        cluster_colors = load_cluster_colors(folder_path)
        if cluster_colors:
            logger.write("Loaded cluster colors from cluster_colors_k6.json\n\n")

        # Find all CSV files with "cluster" in the name
        cluster_files = []
        for file in os.listdir(folder_path):
            if file.endswith('.csv') and 'cluster' in file.lower():
                kmer = get_kmer_from_filename(file)
                if kmer is not None:
                    cluster_files.append((kmer, file))

        if not cluster_files:
            error_msg = "No cluster CSV files found in the specified folder!\n"
            logger.write(error_msg)
            return

        # Sort by k-mer value
        cluster_files.sort(key=lambda x: x[0])

        # Dictionary to store data for each k-mer
        kmer_data = {}
        sequence_to_cluster_map = {}  # Track which cluster each sequence belongs to at each k-mer

        # Process each CSV file
        for kmer, filename in cluster_files:
            file_path = os.path.join(folder_path, filename)

            try:
                df = pd.read_csv(file_path)

                # Validate required columns
                required_cols = ['sequence_id', 'cluster']
                if not all(col in df.columns for col in required_cols):
                    warning_msg = f"Warning: {filename} missing required columns. Skipping...\n"
                    logger.write(warning_msg)
                    continue

                # Count sequences per cluster
                cluster_counts = df['cluster'].value_counts().sort_index()
                total_sequences = len(df)
                unique_clusters = len(cluster_counts)

                # Store the data
                kmer_data[kmer] = {
                    'filename': filename,
                    'cluster_counts': cluster_counts,
                    'total_sequences': total_sequences,
                    'unique_clusters': unique_clusters,
                    'dataframe': df
                }

                # Analyze primer distribution
                primer_stats = analyze_primer_distribution(df)

                # Map sequence IDs to clusters for this k-mer
                sequence_to_cluster_map[kmer] = dict(zip(df['sequence_id'], df['cluster']))

                output = f"K-mer {kmer} ({filename}):\n"
                output += f"  Total sequences: {total_sequences}\n"
                output += f"  Number of clusters: {unique_clusters}\n"
                output += f"  Cluster sizes: {dict(cluster_counts)}\n"

                # Add color information for k6 clusters
                if kmer == 6 and cluster_colors:
                    output += "  Cluster colors:\n"
                    for cluster_id in sorted(cluster_counts.index):
                        cluster_str = str(cluster_id)
                        if cluster_str in cluster_colors:
                            rgba_values = cluster_colors[cluster_str]
                            color_desc = rgba_to_description(rgba_values)
                            output += f"    Cluster {cluster_id}: {color_desc}\n"
                        else:
                            output += f"    Cluster {cluster_id}: No color assigned\n"

                # Add primer distribution per cluster
                output += "  Primer distribution per cluster:\n"
                for cluster_id in sorted(primer_stats.keys()):
                    stats = primer_stats[cluster_id]
                    output += f"    Cluster {cluster_id} ({stats['total_sequences']} sequences):\n"
                    for primer, data in sorted(stats['primers'].items()):
                        output += f"      {primer}: {data['count']} sequences ({data['percentage']:.1f}%)\n"

                output += "\n"
                logger.write(output)

            except Exception as e:
                error_msg = f"Error reading {filename}: {e}\n"
                logger.write(error_msg)
                continue

        # Analyze cluster transitions between k-mer values
        if len(kmer_data) > 1:
            logger.write("\nCluster Transition Analysis:\n")
            logger.write("=" * 50 + "\n")

            kmer_values = sorted(kmer_data.keys())

            for i in range(len(kmer_values) - 1):
                current_kmer = kmer_values[i]
                next_kmer = kmer_values[i + 1]

                logger.write(f"\nTransition from K{current_kmer} to K{next_kmer}:\n")
                logger.write("-" * 30 + "\n")

                current_map = sequence_to_cluster_map[current_kmer]
                next_map = sequence_to_cluster_map[next_kmer]

                # Find common sequences between the two k-mers
                common_sequences = set(current_map.keys()) & set(next_map.keys())

                if not common_sequences:
                    logger.write("No common sequences found between these k-mer values.\n")
                    continue

                # Track how clusters split/merge
                cluster_transitions = defaultdict(lambda: defaultdict(int))

                for seq_id in common_sequences:
                    old_cluster = current_map[seq_id]
                    new_cluster = next_map[seq_id]
                    cluster_transitions[old_cluster][new_cluster] += 1

                summary = f"Common sequences: {len(common_sequences)}\n"
                logger.write(summary)

        # Show cluster evolution tracking at the end
        if len(kmer_data) > 1:
            logger.write("\nCluster Evolution Tracking:\n")
            logger.write("=" * 50 + "\n")

            kmer_values = sorted(kmer_data.keys())

            # Build cluster transition chains
            cluster_chains = {}

            # Start with first k-mer clusters
            first_kmer = kmer_values[0]
            for cluster_id in kmer_data[first_kmer]['cluster_counts'].index:
                cluster_chains[cluster_id] = [cluster_id]

            # For each k-mer transition, find where each cluster's sequences primarily go
            for i in range(len(kmer_values) - 1):
                current_kmer = kmer_values[i]
                next_kmer = kmer_values[i + 1]

                current_map = sequence_to_cluster_map[current_kmer]
                next_map = sequence_to_cluster_map[next_kmer]

                # Find common sequences
                common_sequences = set(current_map.keys()) & set(next_map.keys())

                # For each cluster in current k-mer, find where most of its sequences go
                cluster_transitions = defaultdict(lambda: defaultdict(int))

                for seq_id in common_sequences:
                    old_cluster = current_map[seq_id]
                    new_cluster = next_map[seq_id]
                    cluster_transitions[old_cluster][new_cluster] += 1

                # Update chains based on where majority of sequences go
                new_chains = {}
                for old_cluster in cluster_chains:
                    if old_cluster in cluster_transitions:
                        # Find the destination cluster with the most sequences
                        destinations = cluster_transitions[old_cluster]
                        if destinations:
                            primary_destination = max(destinations, key=destinations.get)
                            sequences_to_primary = destinations[primary_destination]
                            total_sequences = sum(destinations.values())

                            # Extend the chain
                            new_chains[old_cluster] = cluster_chains[old_cluster] + [primary_destination]

                            # Show transition info
                            percentage = (sequences_to_primary / total_sequences) * 100
                            logger.write(
                                f"K{current_kmer} Cluster {old_cluster} → K{next_kmer} Cluster {primary_destination}: {sequences_to_primary}/{total_sequences} sequences ({percentage:.1f}%)\n")

                            # Show splits if significant
                            if len(destinations) > 1:
                                for dest_cluster, count in sorted(destinations.items()):
                                    if dest_cluster != primary_destination and count > 1:
                                        split_percentage = (count / total_sequences) * 100
                                        logger.write(
                                            f"  Also to Cluster {dest_cluster}: {count} sequences ({split_percentage:.1f}%)\n")
                        else:
                            # No sequences found in next k-mer
                            new_chains[old_cluster] = cluster_chains[old_cluster] + ["LOST"]
                    else:
                        # Cluster disappeared
                        new_chains[old_cluster] = cluster_chains[old_cluster] + ["LOST"]

                cluster_chains = new_chains
                logger.write("\n")

            # Show final evolution paths
            logger.write("Final cluster evolution paths:\n")
            for start_cluster in sorted(cluster_chains.keys()):
                path = " > ".join(map(str, cluster_chains[start_cluster]))
                logger.write(f"  Cluster {start_cluster}: {path}\n")

    except Exception as e:
        error_msg = f"\nUnexpected error: {e}\n"
        logger.write(error_msg)

    finally:
        # Program footer with timing
        end_time = datetime.now()
        processing_time = end_time - start_time

        footer = f"""
{'=' * 60}
End Time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}
Total Processing Time: {processing_time}
Script: {script_name}
{'=' * 60}
"""
        logger.write(footer)
        logger.close()

        print(f"\nLog file saved to: {log_file_path}")


if __name__ == "__main__":
    analyze_barcode_clusters()