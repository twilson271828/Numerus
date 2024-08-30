#include "../include/BigInt.hpp"
#include <gtest/gtest.h>

class BigIntTest : public ::testing::Test {
public:
  void SetUp() override {
    //long constructors
  
    //zero_long=BigInt(0);
    //pi_long=BigInt(314159);
    //mpi_long=BigInt(-314159);

    //string constructors
    nan  = BigInt("NaN");
    //c0 = BigInt("0");
    e = BigInt("271828");
    //pi = BigInt("314159");
    me = BigInt("-271828");
    //mpi = BigInt("-314159");
  }
  void TearDown() override {}

  BigInt zero_long;
  BigInt pi_long;
  BigInt mpi_long;

  BigInt nan;
  BigInt c0;
  BigInt e;
  BigInt pi;
  BigInt me;
  BigInt mpi;
};

TEST_F(BigIntTest, ConstructorTests) {

  BigInt z1("NaN");
  BigInt z2("32$$%");

  
  //BigInt z2("314159");
  BigInt z3("271828");
  BigInt z4("-271828");
  //BigInt z4("-314159");
  //BigInt z5("0");

  EXPECT_EQ(z1, nan);
  EXPECT_NE(z2,nan);
  //EXPECT_EQ(z1.get_sign(), UNDEFINED);
  //EXPECT_EQ(z2, pi);
  EXPECT_EQ(z3, e);
  EXPECT_EQ(z4, me);
  //EXPECT_EQ(z4, mpi);
  //EXPECT_EQ(z5, c0);


  //BigInt z6(0);
  //BigInt z7(314159);
  //BigInt z8(-314159);

  //EXPECT_EQ(z6, zero_long);
  //EXPECT_EQ(z7, pi_long);
  //EXPECT_EQ(z8, mpi_long);
}

#if 0
TEST_F (BigIntTest, MultiplicationTests) {
  BigInt z1("271828453454345545545545");
  BigInt z2("314159453453523442343");
  
  BigInt z3 = z1*z2;
  BigInt truth1("85397478370333732998963744994908858288011935");
  EXPECT_EQ(z3, truth1);

  BigInt z4("34232432785792578758923757278528379785274378427843747243823443247224378437892347984378947894378439423789423789423789");
  BigInt z5("234893209423904767780099988797989666766767876767768768");
  BigInt z6=z4*z5;
  BigInt truth2("8040966003442919903678123175346087316928833848925534814654901457889010849276590692709724780423965243729299801473633750753659491030666477046084629010259749637957910421952");

  EXPECT_EQ(z6, truth2);

  BigInt z7("271828453454345545545545");
  BigInt z8("314159453453523442343"); 
  BigInt z9 = z7*z8;
  BigInt truth3("85397478370333732998963744994908858288011935");
  EXPECT_EQ(z9, truth3);

  BigInt z10("-271828453454345545545545");
  BigInt z11("314159453453523442343"); 
  BigInt z12 = z10*z11;
  BigInt truth4("-85397478370333732998963744994908858288011935");
  EXPECT_EQ(z12, truth4);

  BigInt z13("-271828453454345545545545");
  BigInt z14("-314159453453523442343");
  BigInt z15 = z13*z14;
  BigInt truth5("85397478370333732998963744994908858288011935");
  EXPECT_EQ(z15, truth5);



}
TEST_F(BigIntTest, InequalityTests) {

  bool z1 = pi > e;
  EXPECT_EQ(z1, 1);

  bool z2 = e > pi;
  EXPECT_EQ(z2, 0);

  bool z3 = mpi > me;
  EXPECT_EQ(z3, 0);

  bool z4 = me > mpi;
  EXPECT_EQ(z4, 1);

  bool z5 = e > mpi;
  EXPECT_EQ(z5, 1);

  bool z6 = pi < e;
  EXPECT_EQ(z6, 0);

  bool z7 = e < pi;
  EXPECT_EQ(z7, 1);

  bool z8 = mpi < me;
  EXPECT_EQ(z8, 1);

  bool z9 = me < mpi;
  EXPECT_EQ(z9, 0);

  bool z10 = e < pi;
  EXPECT_EQ(z10, 1);

  bool z11 = e < mpi;
  EXPECT_EQ(z11, 0);
}

TEST_F(BigIntTest, EqualityTests) { 
  EXPECT_EQ(e, e); 
  EXPECT_NE(e, pi);
  }

// Demonstrate some basic assertions.
TEST_F(BigIntTest, AdditionTests) {

  BigInt z1("4535435723");
  BigInt z2("34355");
  BigInt mz1("-4535435723");
  BigInt mz2("-34355");

  // x > 0, y > 0
  BigInt truth1("4535470078");
  EXPECT_EQ(z1 + z2, truth1);

  // x < 0, y < 0
  BigInt truth2("-4535470078");
  EXPECT_EQ(mz1 + mz2, truth2);

  // x > 0, y < 0
  BigInt truth3("4535401368");
  EXPECT_EQ(z1 + mz2, truth3);

  // x < 0, y > 0
  BigInt truth4("-4535401368");
  EXPECT_EQ(mz1 + z2, truth4);

  BigInt truth6("0");
  EXPECT_EQ(z1 + mz1, truth6);
}

TEST_F(BigIntTest, shift_n_Test) {

  BigInt z("271828");
  BigInt z3a("000271828");
  BigInt z3b("271828000");
  BigInt z_add_3_to_front = z.shift_n(3, true);
  BigInt z_add_3_to_rear = z.shift_n(3);
  EXPECT_EQ(z3a, z_add_3_to_front);
  EXPECT_EQ(z3b, z_add_3_to_rear);
}

TEST_F(BigIntTest, split_it_tests) {
  
  
  BigInt z("6789424643665457123213125523442134324234242352342380724234242");
  split c1 = z.split_it(6);
  BigInt zleft("6789424643665457123213125523442134324234242352342380724");
  BigInt zright("234242");
  BigInt nan("NaN");
  EXPECT_EQ(c1.xleft, zleft);
  EXPECT_EQ(c1.xright, zright);
  split c2 = z.split_it(62);
  EXPECT_EQ(c2.xleft,nan);
  EXPECT_EQ(c2.xright,nan);

}

TEST_F(BigIntTest, SliceTests) {

  BigInt z("6789353555355553535353553535");
  BigInt zslice1("67893535553");
  BigInt z1 = z.slice(0, 10);
  EXPECT_EQ(zslice1, z1);
  BigInt z2 = z.slice(-1, -1);
  EXPECT_EQ(z2.get_sign(), UNDEFINED);

  BigInt z3 = z.slice(5, 5);
  EXPECT_EQ(z3, BigInt("5"));

  BigInt z4 = z.slice(0, 1000);
  EXPECT_EQ(z4.get_sign(), UNDEFINED);

  
  BigInt z5 = z.slice(27, 34);
  EXPECT_EQ(z5.get_sign(), UNDEFINED);

  BigInt z6 = z.slice(27, 28);
  EXPECT_EQ(z6.get_sign(), UNDEFINED);

  BigInt z7 = z.slice(z.size() - 1, z.size() - 1);
  EXPECT_EQ(z7, BigInt("5"));
  
}

TEST_F(BigIntTest, SubtractionTests) {

  BigInt z1("4324432423432");
  BigInt mz1("-4324432423432");

  BigInt z2("23245");
  BigInt mz2("-23245");

  // x > y, x > 0, y>0
  BigInt truth1("4324432400187");
  EXPECT_EQ(z1 - z2, truth1);

  // x > y, x> 0,y < 0
  BigInt truth2("4324432446677");
  EXPECT_EQ(z1 - mz2, truth2);

  // x > y, x < 0,y < 0
  BigInt truth3("4324432400187");
  EXPECT_EQ(mz2 - mz1, truth3);

  // x < y, x < 0, y < 0
  BigInt truth4("-4324432400187");
  EXPECT_EQ(mz1 - mz2, truth4);

  // x < y, x > 0,y > 0
  BigInt truth5("-4324432400187");
  EXPECT_EQ(z2 - z1, truth5);

  BigInt truth6("0");
  EXPECT_EQ(z2 - z2, truth6);
  // EXPECT_EQ(truth4,mz1-mz2);
}
#endif

int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
