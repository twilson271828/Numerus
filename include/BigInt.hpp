    #include <iostream>
    #include <cstdlib>
    #include <vector>
    #include <string>
    #include <cmath>
    #include <complex>
    //#include <armadillo>
    
    enum SIGN {POS,NEG};

    class BigInt {

        private:

         std::vector<int> numerus;

         SIGN sign;

        private:
            BigInt vsub(BigInt &x,BigInt &y) const;
            BigInt vadd(BigInt &x,BigInt &y) const;
            BigInt vmult(BigInt &x, BigInt &y) const;
            BigInt karatsuba(BigInt &x, BigInt &y) const;
            BigInt Schonhage_Strassen(BigInt &x,BigInt&y) const;
            BigInt Toom3(BigInt &x, BigInt &y) const;
            std::complex<double> exponentiate(size_t k,size_t n,size_t N);
            std::complex<double> dift(std::vector<std::complex<double>> &input,size_t n); 
            std::complex<double> dft(std::vector<std::complex<double>> &input,size_t n); 
            size_t bitrev(size_t n);
        public: 

        BigInt();
        
        /// @brief 
        /// @param c 
        BigInt(const std::string c);

        BigInt(const BigInt &num);

        BigInt(const long &num);

        /// @brief 
        /// @param m 
        /// @param add_to_front 
        BigInt m10(const int m, bool add_to_front = false) const;
       
        int operator[](const int i) const;
        int size() const;

        int get_sign() const;

        /// @brief 
        void negative();

        
        void insert(const int &val,const int &ix);
        
        
        BigInt operator *(const BigInt &num);

        /// @brief 
        /// @param num 
        /// @return 
        
        BigInt operator /(const BigInt &num) const;

        BigInt operator+(const BigInt &num) const;

        BigInt operator -(const BigInt &num) const;

        bool operator < (const BigInt &num) const;

        bool operator <= (const BigInt &num) const;

        bool operator >= (const BigInt &num) const;

        bool operator > (const BigInt &num) const;

        bool operator != (const BigInt &num) const;

        bool operator == (const BigInt &num) const;

        BigInt operator ++ ();

        BigInt operator --();

        BigInt operator ! ()const; 


        
        friend std::ostream& operator<<(std::ostream &out,const BigInt& num);

        

    };