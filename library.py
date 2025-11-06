from fractions import Fraction
from typing import Optional


def convert_to_triangular_form(matrix: list) -> list:
    """
    Attempts to transform a matrix into a triangular form.
    Resource: https://en.wikipedia.org/wiki/Triangular_matrix
    :param matrix: The matrix must be a valid Augmented Matrix in the form n x m matrix.
    :return: The triangular form of the matrix.
    :throws: ValueError if the matrix does not meet the requirements for an augmented matrix (n x m, consisting of an n x n coefficient matrix and an n x 1 constant vector).
    :throws: ArithmeticError if the matrix cannot be transformed into a triangular form.
    """
    # break the matrix into its coefficient matrix and constant vector to validate shape
    coefficient_matrix = [row[:-1] for row in matrix]
    constant_vector = [row[-1] for row in matrix]

    # Check if the coefficient matrix is n x n
    if len(coefficient_matrix) != len(coefficient_matrix[0]):
        raise ValueError("The matrix must be square.")
    # Check if the constant vector is n x 1
    if len(constant_vector) != len(coefficient_matrix):
        raise ValueError(f"The vector must be a vector of length {len(coefficient_matrix)}.")

    # Initialize the triangular form.
    aM = [[*coefficient_matrix[i], constant_vector[i]] for i in range(len(coefficient_matrix))]

    # Iterate through the rows.
    for t in range(len(aM) -1):
        # Step 1: Rearrange rows t through len(aM) - 1 so that the element in column t is not zero, else raise ArithemeticError
        s = 0
        while aM[t][t] == 0 and s < len(aM) - t:
            s += 1
            # Move the row to the bottom of the matrix.
            aM.append(aM.pop(t))

        if s == len(aM) - t:
            raise ArithmeticError("The matrix cannot be transformed into a triangular form.")

        # Step 2: for rows t + 1 through len(aM) - 1, Add a nonZero multiple of row t to row t + 1 inorder to make the element in column t zero.
        for i in range(t + 1, len(aM)):
            if aM[i][t] != 0:
                # Calculate the nonZero multiple of row t.
                factor = aM[i][t] / aM[t][t]
                # Calculate the adjusted row t that will be used to make the element in column t zero.
                row = [x * factor for x in aM[t]]
                # replace the first t elements of the row with zeroes so that we are only adding to the last t elements.
                # row = [0] * t + row[t:]
                # Replace aM[i][t] with aM[i][t] + row
                aM[i] = [y - x for x, y in zip(row, aM[i])]

    # Return the triangular form.
    return aM

def find_solution_with_triangular_form(matrix: list) -> Optional[list[int]]:
    """
    Given a system of equations as an augmented matrix, n x n + n x r, return the solution as a vector.
    :param matrix: An m x n matrix representing a system of equations
    :return: A vector representing the solution to the system of equations.
    """
    try:
        tM = convert_to_triangular_form(matrix)
    except ArithmeticError:
        return None
    except ValueError:
        return None

    # Step 2: Apply Back Substitution to the triangular form.
    solution = [0] * len(tM)
    for i in range(len(tM) - 1, -1, -1):
        # Create a row by substituting the solution for each previous term.
        r = []
        for j in range(len(tM) + 1):
            if j < i:
                r.append(0)
            elif j == i:
                r.append(tM[i][j])
            elif j == len(tM):
                r.append(tM[i][j])
            else:
                r.append(tM[i][j] * solution[j])

        # Now Subtract each previous term from the constant term.
        for j in range(i + 1, len(tM)):
            r[len(tM)] -= r[j]

        # Now divide the constant term by the coefficient of the i term.
        solution[i] = r[len(tM)] / tM[i][i]
    return solution

def to_mixed_fraction_string(f: Fraction) -> str:
    """
    Converts a Fraction to a mixed number string (e.g., '1 1/3' or '2').
    :param f: The Fraction to convert.
    :return: The string representation of the Fraction.
    """
    if f.denominator == 1:
        return str(f.numerator)

    whole = f.numerator // f.denominator
    remainder = f.numerator % f.denominator

    if whole == 0 and remainder != 0:
        return f"{remainder}/{f.denominator}"
    elif remainder == 0:
        return str(whole)
    else:
        # Handle negative mixed numbers properly (e.g. -1 1/3)
        if whole < 0:
            return f"{whole} {abs(remainder)}/{f.denominator}"
        return f"{whole} {remainder}/{f.denominator}"

def find_and_format_solution_to_a_system_of_equations_via_conversion_to_triangular_form(matrix: list) -> str:
    """
    Solves a system of equations using the triangular form method and formats the solution.
    :param matrix: An m x n matrix representing a system of equations
    :return: The formatted solution for the system of equations.
    """
    if solution := find_solution_with_triangular_form(matrix):
        # Convert each number to a string, handling integer floats
        formatted_solution = []
        for num in solution:
            # Convert the float to a Fraction, limiting denominator complexity
            # to avoid overly long fractions due to float precision issues.
            # A max_denominator of 100 often gives human-readable results.
            frac = Fraction(num).limit_denominator(100)
            formatted_solution.append(to_mixed_fraction_string(frac))

        # Join the list of strings
        return f"solution = ({', '.join(formatted_solution)})"
    return "No solution."

def solve_with_triangular_form(matrix: list) -> str:
    """
    An alias for find_and_format_solution_to_a_system_of_equations_via_conversion_to_triangular_form.
    :param matrix: An m x n matrix representing a system of equations
    :return: The formatted solution for the system of equations.
    """
    return find_and_format_solution_to_a_system_of_equations_via_conversion_to_triangular_form(matrix)
