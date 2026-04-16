# This script extracts the number of correct coefficients from a specific line in a text file.
import sys

def extract_correct_coefficients(file_path, line_number):
    """
    Extracts the number of correct coefficients from the specified line in the file.

    Args:
        file_path (str): Path to the text file.
        line_number (int): Line number to extract (1-based index).

    Returns:
        int: The extracted number of correct coefficients.
    """
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            if line_number > len(lines) or line_number <= 0:
                raise ValueError("Line number is out of range.")

            # Extract the specific line
            target_line = lines[line_number - 1].strip()

            # Extract the number before "correct coefficients"
            if "correct coefficients" in target_line:
                #parts = target_line.split("recovered coefficients")
                parts = target_line.split("correct coefficients")
                number = parts[0].strip().split()[-1]
                return int(number)
            else:
                raise ValueError("The line does not contain 'correct coefficients'.")

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
    except ValueError as ve:
        print(f"Error: {ve}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Example usage
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python extract_correct_coefficients.py <directory>")
        sys.exit(1)

    directory = sys.argv[1]

    for p in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        file_path = f"{directory}/{p}.txt"  # Use the directory argument
        line_number = 371  # Replace with the desired line number
        result = extract_correct_coefficients(file_path, line_number)
        if result is not None:
            print(f"Extracted number of correct coefficients for p={p}: {result}")