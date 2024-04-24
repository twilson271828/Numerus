#include "../include/BigInt.hpp"
#include <gtest/gtest.h>

class BigIntTest : public ::testing::Test {
public:
  void SetUp() override {

    c0 = BigInt("0");
    e = BigInt("271828");
    pi = BigInt("314159");
    me = BigInt("-271828");
    mpi = BigInt("-314159");
  }
  void TearDown() override {}

  BigInt c0;
  BigInt e;
  BigInt pi;
  BigInt me;
  BigInt mpi;
};


TEST_F(BigIntTest,InequalityTests) {

   
 bool z1 = pi > e;
 EXPECT_EQ(z1,1);

 bool z2 = e > pi;
 EXPECT_EQ(z2,0);
 
 bool z3 = mpi > me;
 EXPECT_EQ(z3,0);

bool z4 = me > mpi;
EXPECT_EQ(z4,1) ;

bool z5 = e > mpi;
EXPECT_EQ(z5,1);

bool z6 = pi < e;
EXPECT_EQ(z6,0);

bool z7 = e < pi;
EXPECT_EQ(z7,1);

bool z8 = mpi < me;
EXPECT_EQ(z8,1);

bool z9 = me < mpi;
EXPECT_EQ(z9, 0);

bool z10 = e < pi;
EXPECT_EQ(z10,1);

bool z11 = e < mpi;
EXPECT_EQ(z11, 0);

}


TEST_F(BigIntTest,EqualityTests) {

  EXPECT_EQ(e,e);

}


// Demonstrate some basic assertions.
TEST_F(BigIntTest, AdditionTests) {

  BigInt z1("4535435723");
  BigInt z2("34355");
  BigInt mz1("-4535435723");
  BigInt mz2("-34355");

  //x > 0, y > 0
  BigInt truth1("4535470078");
  EXPECT_EQ(z1+z2,truth1);

  //x < 0, y < 0
  BigInt truth2("-4535470078");
  EXPECT_EQ(mz1+mz2,truth2);

  //x > 0, y < 0
  BigInt truth3("4535401368");
  EXPECT_EQ(z1+mz2,truth3);

   //x < 0, y > 0
  BigInt truth4("-4535401368");
  EXPECT_EQ(mz1+z2,truth4);

  BigInt truth6("0");
  EXPECT_EQ(z1+mz1,truth6);
  
}

TEST_F(BigIntTest,shift_n_Test) {

  BigInt z("271828");
  BigInt z3a("000271828");
  BigInt z3b("271828000");
  BigInt z_add_3_to_front = z.shift_n(3,true);
  BigInt z_add_3_to_rear = z.shift_n(3);
  EXPECT_EQ(z3a,z_add_3_to_front);
  EXPECT_EQ(z3b,z_add_3_to_rear);

}

TEST_F(BigIntTest,split_it_tests) {

BigInt z("6789424643665457123213125523442134324234242352342380724234242");
EXPECT_EQ(z,z);


}

TEST_F(BigIntTest,SliceTests) {

  BigInt z("6789353555355553535353553535");
  BigInt zslice1("67893535553");
  BigInt z1 = z.slice(0,10);
  EXPECT_EQ(zslice1,z1);

  BigInt z2 = z.slice(-1,-1);
  EXPECT_EQ(z2.get_sign(),0);
  BigInt z3 = z.slice(5,5);
  EXPECT_EQ(z3,BigInt("5"));
  BigInt z4 = z.slice(0,1000);
  EXPECT_EQ(z4.get_sign(),0);

  BigInt z5 = z.slice(27,34);
  EXPECT_EQ(z5.get_sign(),0);

  BigInt z6 = z.slice(27,28);
  EXPECT_EQ(z6.get_sign(),0);

  BigInt z7 = z.slice(z.size()-1,z.size()-1);
  EXPECT_EQ(z7,BigInt("5") );
 
}

TEST_F(BigIntTest, SubtractionTests) {
  
  BigInt z1("4324432423432");
  BigInt mz1("-4324432423432");
  
  BigInt z2("23245");
  BigInt mz2("-23245");
  

  //x > y, x > 0, y>0
  BigInt truth1("4324432400187");
  EXPECT_EQ(z1-z2,truth1);

  //x > y, x> 0,y < 0
  BigInt truth2("4324432446677");
  EXPECT_EQ(z1-mz2,truth2);

  //x > y, x < 0,y < 0
  BigInt truth3("4324432400187");
  EXPECT_EQ(mz2 - mz1,truth3);

  //x < y, x < 0, y < 0
   BigInt truth4("-4324432400187");
   EXPECT_EQ(mz1-mz2,truth4);

   //x < y, x > 0,y > 0
   BigInt truth5("-4324432400187");
   EXPECT_EQ(z2-z1,truth5);

  BigInt truth6("0");
  EXPECT_EQ(z2-z2,truth6);
  //EXPECT_EQ(truth4,mz1-mz2);

}


int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
