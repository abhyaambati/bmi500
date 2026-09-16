import random


def dot_product(a, b):
    """Compute the dot product of two vectors using a for loop.

    Args:
        a: First vector (sequence of numbers).
        b: Second vector (sequence of numbers), same length as a.

    Returns:
        Scalar sum of element-wise products.

    Raises:
        TypeError: If a or b is not a sequence, or if any element is non-numeric.
        ValueError: If a and b have different lengths, or either is empty.
    """
    # Guard: inputs must be sequences (list, tuple, etc.), not scalars/strings alone
    if not hasattr(a, "__len__") or not hasattr(b, "__len__"):
        raise TypeError("Both inputs must be sequences (e.g. lists or tuples)")
    if isinstance(a, (str, bytes)) or isinstance(b, (str, bytes)):
        raise TypeError("Inputs must be numeric sequences, not strings")

    # Guard: empty vectors have no well-defined non-trivial product for this assignment
    if len(a) == 0 or len(b) == 0:
        raise ValueError("Vectors must be non-empty")

    # Guard: lengths must match for a valid dot product
    if len(a) != len(b):
        raise ValueError(
            f"Vector length mismatch: len(a)={len(a)}, len(b)={len(b)}"
        )

    result = 0
    # Accumulate a[i] * b[i] for each index using an explicit for loop
    for i in range(len(a)):
        # Guard: reject non-numeric elements early with a clear error
        try:
            result += a[i] * b[i]
        except TypeError as exc:
            raise TypeError(
                f"Non-numeric element at index {i}: a[{i}]={a[i]!r}, b[{i}]={b[i]!r}"
            ) from exc
    return result


def matvec_multiply(matrix, vector):
    """Compute the matrix-vector product using dot_product.

    For an m x n matrix A and length-n vector x, returns the length-m vector
    y where y[i] = A[i] · x.

    Args:
        matrix: 2D sequence of numbers (list of rows).
        vector: 1D sequence of numbers whose length equals the number of columns.

    Returns:
        List of floats/ints: the matrix-vector product.

    Raises:
        TypeError: If matrix/vector structure or element types are invalid.
        ValueError: If dimensions are incompatible or inputs are empty.
    """
    # Guard: matrix must be a non-empty sequence of rows
    if not hasattr(matrix, "__len__") or isinstance(matrix, (str, bytes)):
        raise TypeError("matrix must be a sequence of row sequences")
    if len(matrix) == 0:
        raise ValueError("matrix must be non-empty")

    # Guard: vector must be a non-empty numeric sequence
    if not hasattr(vector, "__len__") or isinstance(vector, (str, bytes)):
        raise TypeError("vector must be a numeric sequence")
    if len(vector) == 0:
        raise ValueError("vector must be non-empty")

    n_cols = len(vector)
    result = []
    # Each output entry is the dot product of one matrix row with the vector
    for row_idx, row in enumerate(matrix):
        if not hasattr(row, "__len__") or isinstance(row, (str, bytes)):
            raise TypeError(f"matrix row {row_idx} must be a numeric sequence")
        # Guard: every row must have the same length as the vector (column count)
        if len(row) != n_cols:
            raise ValueError(
                f"Row {row_idx} has length {len(row)}, expected {n_cols} "
                f"(matching vector length)"
            )
        result.append(dot_product(row, vector))
    return result


def main():
    """Test matvec_multiply with randomly generated 1000x1000 data.

    Builds a 1000x1000 matrix and a length-1000 vector with uniform random
    floats in [0, 1), runs the product, and prints a short summary so the
    script can be run as a smoke test: python matvec_multiply.py
    """
    n = 1000
    # Random matrix A (n x n) and vector x (n,)
    matrix = [[random.random() for _ in range(n)] for _ in range(n)]
    vector = [random.random() for _ in range(n)]

    product = matvec_multiply(matrix, vector)

    print(f"Matrix shape: {n} x {n}")
    print(f"Vector length: {n}")
    print(f"Result length: {len(product)}")
    print(f"Result[0] (sample): {product[0]:.6f}")
    print(f"Result sum (sanity check): {sum(product):.6f}")


if __name__ == "__main__":
    main()
