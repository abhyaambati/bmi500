"""Tests for matrix-vector and dot-product functions in matvec_multiply.py."""

import unittest

from matvec_multiply import dot_product, matvec_multiply


class TestDotProduct(unittest.TestCase):
    """Tests for the dot_product function."""

    def test_basic_integers(self):
        self.assertEqual(dot_product([1, 2, 3], [4, 5, 6]), 32)

    def test_floats(self):
        self.assertAlmostEqual(dot_product([1.5, 2.0], [2.0, 0.5]), 4.0)

    def test_single_element(self):
        self.assertEqual(dot_product([7], [3]), 21)

    def test_zeros(self):
        self.assertEqual(dot_product([0, 0, 0], [1, 2, 3]), 0)

    def test_tuples_accepted(self):
        self.assertEqual(dot_product((1, 2), (3, 4)), 11)

    def test_length_mismatch_raises(self):
        with self.assertRaises(ValueError):
            dot_product([1, 2], [1, 2, 3])

    def test_empty_raises(self):
        with self.assertRaises(ValueError):
            dot_product([], [])

    def test_non_numeric_raises(self):
        with self.assertRaises(TypeError):
            dot_product([1, "a"], [2, 3])

    def test_string_input_raises(self):
        with self.assertRaises(TypeError):
            dot_product("ab", "cd")


class TestMatvecMultiply(unittest.TestCase):
    """Tests for the matvec_multiply function."""

    def test_identity(self):
        identity = [[1, 0], [0, 1]]
        vector = [3, 4]
        self.assertEqual(matvec_multiply(identity, vector), [3, 4])

    def test_known_product(self):
        matrix = [[1, 2, 3], [4, 5, 6]]
        vector = [1, 0, -1]
        # row0: 1*1 + 2*0 + 3*(-1) = -2
        # row1: 4*1 + 5*0 + 6*(-1) = -2
        self.assertEqual(matvec_multiply(matrix, vector), [-2, -2])

    def test_single_row(self):
        self.assertEqual(matvec_multiply([[2, 3]], [4, 5]), [23])

    def test_column_vector_style(self):
        # 3x1 matrix times length-1 vector
        matrix = [[2], [3], [4]]
        self.assertEqual(matvec_multiply(matrix, [10]), [20, 30, 40])

    def test_row_length_mismatch_raises(self):
        matrix = [[1, 2], [3]]  # jagged
        with self.assertRaises(ValueError):
            matvec_multiply(matrix, [1, 2])

    def test_dimension_mismatch_raises(self):
        matrix = [[1, 2, 3], [4, 5, 6]]
        with self.assertRaises(ValueError):
            matvec_multiply(matrix, [1, 2])  # vector too short

    def test_empty_matrix_raises(self):
        with self.assertRaises(ValueError):
            matvec_multiply([], [1, 2])

    def test_empty_vector_raises(self):
        with self.assertRaises(ValueError):
            matvec_multiply([[1, 2]], [])

    def test_non_numeric_matrix_raises(self):
        with self.assertRaises(TypeError):
            matvec_multiply([[1, "x"], [2, 3]], [1, 1])

    def test_invalid_matrix_type_raises(self):
        with self.assertRaises(TypeError):
            matvec_multiply("not a matrix", [1, 2])


if __name__ == "__main__":
    unittest.main()
