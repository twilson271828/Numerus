import unittest
import random
from gemini_schonhage_strassen import SchonhageStrassen
# from your_script_name import SchonhageStrassen  # Use this if in a separate file

class TestSchonhageStrassen(unittest.TestCase):
    def setUp(self):
        """Initialize the algorithm instance before each test."""
        self.ssa = SchonhageStrassen()

    def test_basic_multiplication(self):
        """Test simple multiplication to ensure baseline functionality."""
        a = 123
        b = 456
        expected = a * b
        result = self.ssa.multiply(a, b)
        self.assertEqual(result, expected, f"Failed basic: {a} * {b}")

    def test_different_lengths_small_large(self):
        """
        Test multiplying a very small number by a very large number.
        This verifies that the polynomial padding (A += [0]...) works correctly.
        """
        small = 15
        # 2 to the power of 100 (approx 31 digits)
        large = 1 << 100 
        
        expected = small * large
        result = self.ssa.multiply(small, large)
        
        self.assertEqual(result, expected, "Failed multiplying small * large integer")

    def test_different_lengths_massive_difference(self):
        """
        Test extremely disparate lengths (e.g., 1 digit vs 1000 digits).
        """
        tiny = 7
        # Create a massive number: 1 followed by 300 zeros
        huge = int('1' + '0' * 300) 
        
        expected = tiny * huge
        result = self.ssa.multiply(tiny, huge)
        self.assertEqual(result, expected, "Failed massive length difference test")

    def test_large_equal_length(self):
        """Test multiplying two large numbers of roughly equal length."""
        # Two 50-digit numbers
        a = 83748237489237489237489237489237489237423894723894
        b = 93847239847239847239847239847298374982374982374982
        
        expected = a * b
        result = self.ssa.multiply(a, b)
        self.assertEqual(result, expected, "Failed large equal length test")

    def test_zero_multiplication(self):
        """Test multiplication by zero (edge case)."""
        large = 1234567890123456789
        self.assertEqual(self.ssa.multiply(large, 0), 0)
        self.assertEqual(self.ssa.multiply(0, large), 0)

    def test_identity_multiplication(self):
        """Test multiplication by 1."""
        val = 987654321
        self.assertEqual(self.ssa.multiply(val, 1), val)
        self.assertEqual(self.ssa.multiply(1, val), val)

    def test_negative_numbers(self):
        """Test logic for handling signs."""
        a = -12345
        b = 67890
        
        # Negative * Positive
        self.assertEqual(self.ssa.multiply(a, b), a * b)
        # Negative * Negative
        self.assertEqual(self.ssa.multiply(a, -b), a * (-b))

    def test_fuzzy_random_lengths(self):
        """
        Generate random pairs of integers with random lengths to find edge cases.
        This runs 100 randomized tests.
        """
        print("\nRunning fuzzy tests (100 iterations)...")
        for _ in range(100):
            # Random bit lengths between 1 and 2048 bits
            bits_a = random.randint(1, 2048)
            bits_b = random.randint(1, 2048)
            
            a = random.getrandbits(bits_a)
            b = random.getrandbits(bits_b)
            
            # Occasionally make one number negative
            if random.random() < 0.2: a = -a
            if random.random() < 0.2: b = -b
            
            expected = a * b
            result = self.ssa.multiply(a, b)
            
            self.assertEqual(result, expected, 
                             f"Fuzzy test failed for inputs of length {bits_a} and {bits_b} bits")
        print("Fuzzy tests passed.")

if __name__ == '__main__':
    unittest.main()