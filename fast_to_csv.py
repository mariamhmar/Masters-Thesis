# Import the necessary library
import sys

# Define a function to parse the FASTA file and write to CSV
def fasta_to_csv(input_file, output_file):
    try:
        with open(input_file, 'r') as fasta_file, open(output_file, 'w') as csv_file:
            current_id = None
            sequence = []

            for line in fasta_file:
                line = line.strip()
                if line.startswith('>'):
                    # Write the previous sequence to CSV (if any)
                    if current_id:
                        csv_file.write(f"{current_id},{''.join(sequence)}\n")
                    current_id = line[1:]  # Remove the '>' character
                    sequence = []
                else:
                    sequence.append(line)

            # Write the last sequence to CSV
            if current_id:
                csv_file.write(f"{current_id},{''.join(sequence)}\n")

        print(f"Conversion complete. Output saved to {output_file}")

    except FileNotFoundError:
        print("Input FASTA file not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

# Check if the script is run from the command line
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python fasta_to_csv.py input.fasta output.csv")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        fasta_to_csv(input_file, output_file)
