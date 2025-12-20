import numpy as np
import matplotlib.pyplot as plt
from gemini_schonhage_strassen import SchonhageStrassen
import random

def main():
    ssa = SchonhageStrassen()
    small = 1554353453453453535345
    # 2 to the power of 200 (approx 31 digits)
    large = 1 << 200
    
    expected = small * large
    result = ssa.multiply(small, large)

    print(result)
    print(expected)
        




if __name__=="__main__":
    main()
