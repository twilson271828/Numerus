#include "../include/BigInt.hpp"
#include <gtest/gtest.h>

class BigIntTest : public ::testing::Test {
public:
  void SetUp() override {
    // long constructors

    zero_long = BigInt(0);
    pi_long = BigInt(314159);
    mpi_long = BigInt(-314159);

    // string constructors
    nan1 = BigInt("NaN");
    nan2 = BigInt("32$$%");
    nan3 = BigInt("%");

    c0 = BigInt("0");
    e = BigInt("271828");
    pi = BigInt("314159");
    me = BigInt("-271828");
    mpi = BigInt("-314159");
  }
  void TearDown() override {}

  BigInt zero_long;
  BigInt pi_long;
  BigInt mpi_long;
  BigInt empty;
  BigInt nan1;
  BigInt nan2;
  BigInt nan3;
  BigInt c0;
  BigInt e;
  BigInt pi;
  BigInt me;
  BigInt mpi;
};

TEST_F(BigIntTest, LongConstructorTests) {
  EXPECT_EQ(zero_long, c0);
  EXPECT_EQ(pi_long, pi);
  EXPECT_EQ(mpi_long, mpi);
}

TEST_F(BigIntTest, print_numerus) {
  BigInt z("271828");
  z.print_numerus();
}

TEST_F(BigIntTest, ConstructorTests) {

  BigInt z1("NaN");
  BigInt z2("32$$%");

  BigInt z3("271828");
  BigInt z4("-271828");

  BigInt z5("0");
  BigInt z6("+314159");
  BigInt z7("+-314159");
  BigInt z8("-+314159");
  BigInt z9("&");
  BigInt z10("&&3434324");

  EXPECT_EQ(z1, nan1);
  EXPECT_NE(z2, nan1);
  EXPECT_EQ(z1.get_sign(), UNDEFINED);
  EXPECT_EQ(z9, nan1);
  EXPECT_EQ(z9.get_sign(), UNDEFINED);

  EXPECT_EQ(z10, nan1);
  EXPECT_EQ(z10.get_sign(), UNDEFINED);

  EXPECT_EQ(z3, e);
  EXPECT_EQ(z4, me);
  EXPECT_EQ(z5, c0);
  EXPECT_EQ(z6, pi);
  EXPECT_EQ(z7, nan1);
  EXPECT_EQ(z8, nan1);
}

TEST_F(BigIntTest, Get_Numerus) {
  BigInt z1("271828");
  std::vector<uint8_t> v1 = {2, 7, 1, 8, 2, 8};
  EXPECT_EQ(z1.getNumerus(), v1);
}

TEST_F(BigIntTest, OstreamOperator) {

  std::ostringstream oss;
  BigInt z;
  std::vector<uint8_t> v1 = {0, 0, 0, 2, 7, 1, 8, 2, 8};
  z.setNumerus(v1);
  z.set_sign(POS);
  oss << z;
  EXPECT_EQ(oss.str(), "271828");

  oss.str(""); // Clear the stream
  BigInt zempty;
  std::vector<uint8_t> vempty;
  zempty.setNumerus(vempty);
  oss << zempty;
  EXPECT_EQ(oss.str(), "_NULL");

  oss.str(""); // Clear the stream
  oss << e;
  EXPECT_EQ(oss.str(), "271828");

  oss.str(""); // Clear the stream
  oss << pi;
  EXPECT_EQ(oss.str(), "314159");

  oss.str(""); // Clear the stream
  oss << me;
  EXPECT_EQ(oss.str(), "-271828");

  oss.str(""); // Clear the stream
  oss << mpi;
  EXPECT_EQ(oss.str(), "-314159");

  oss.str("");
  oss << nan1;
  EXPECT_EQ(oss.str(), "NaN");

  oss.str("");
  oss << nan2;
  EXPECT_EQ(oss.str(), "NaN");

  oss.str("");
  oss << empty;
  EXPECT_EQ(oss.str(), "_NULL");
}

TEST_F(BigIntTest, MultiplicationTests) {
  BigInt z1("271828453454345545545545");
  BigInt z2("314159453453523442343");

  BigInt z3 = z1 * z2;

  BigInt truth1("85397478370333732998963744994908858288011935");
  EXPECT_EQ(z3, truth1);

  BigInt z4("342324327857925787589237572785283797852743784278437472438234432472"
            "24378437892347984378947894378439423789423789423789");
  BigInt z5("234893209423904767780099988797989666766767876767768768");
  BigInt z6 = z4 * z5;
  BigInt truth2("80409660034429199036781231753460873169288338489255348146549014"
                "57889010849276590692709724780423965243729299801473633750753659"
                "491030666477046084629010259749637957910421952");

  EXPECT_EQ(z6, truth2);

  BigInt z7("271828453454345545545545");
  BigInt z8("314159453453523442343");
  BigInt z9 = z7 * z8;
  BigInt truth3("85397478370333732998963744994908858288011935");
  EXPECT_EQ(z9, truth3);

  BigInt z10("-271828453454345545545545");
  BigInt z11("314159453453523442343");
  BigInt z12 = z10 * z11;
  BigInt truth4("-85397478370333732998963744994908858288011935");
  EXPECT_EQ(z12, truth4);

  BigInt z13("-271828453454345545545545");
  BigInt z14("-314159453453523442343");
  BigInt z15 = z13 * z14;
  BigInt truth5("85397478370333732998963744994908858288011935");
  EXPECT_EQ(z15, truth5);

  BigInt z16("535344");
  BigInt z17("42323");
  BigInt truth6("22657364112");
  BigInt z18 = z16 * z17;
  EXPECT_EQ(z18, truth6);

  BigInt z19;
  BigInt z20;
  BigInt z21 = z19 * z20;
  EXPECT_EQ(z21, nan1);

  BigInt z22("0");
  BigInt z23("0");
  BigInt z24 = z22 * z23;
  EXPECT_EQ(z24, c0);

  BigInt z25(0);
  BigInt z26(5);
  BigInt z27 = z25 * z26;
  EXPECT_EQ(z27, c0);

  BigInt z28(34242);
  BigInt z29(432546366);
  BigInt z30 = z28 * z29;
  EXPECT_EQ(z30, BigInt("14811252664572"));
}

TEST_F(BigIntTest, DivisionTests) {
  BigInt x("7294372378472835723758");
  BigInt y("2568");

  divmod10 z = x / y;
  BigInt truth_quotient("2840487686321197711");
  BigInt truth_remainder("1910");

  EXPECT_EQ(z.quotient, truth_quotient);
  EXPECT_EQ(z.remainder, truth_remainder);
}

TEST_F(BigIntTest, LEQ) {
  // n == m
  BigInt e("271828");
  BigInt pi("314159");

  BigInt me("-271828");
  BigInt mpi("-314159");

  bool t1 = e <= pi;
  EXPECT_EQ(t1, 1);

  bool t2 = pi <= e;
  EXPECT_EQ(t2, 0);

  bool t3 = me <= mpi;
  EXPECT_EQ(t3, 0);

  bool t4 = mpi <= me;
  EXPECT_EQ(t4, 1);

  bool t5 = e <= e;
  EXPECT_EQ(t5, 1);

  bool t6 = e <= me;
  EXPECT_EQ(t6, 0);

  bool t7 = me <= e;
  EXPECT_EQ(t7, 1);

  // n > m
  BigInt pipi("-314159654");

  bool t8 = pipi <= mpi;
  EXPECT_EQ(t8, 1);

  bool t9 = pi <= e;
  EXPECT_EQ(t9, 0);

  bool t10 = e <= pi;
  EXPECT_EQ(t10, 1);
}

TEST_F(BigIntTest, GEQ) {
  // n == m
  BigInt e("271828");
  BigInt pi("314159");

  BigInt me("-271828");
  BigInt mpi("-314159");

  bool t1 = e >= pi;
  EXPECT_EQ(t1, 0);

  bool t2 = pi >= e;
  EXPECT_EQ(t2, 1);

  bool t3 = me >= mpi;
  EXPECT_EQ(t3, 1);

  bool t4 = mpi >= me;
  EXPECT_EQ(t4, 0);

  bool t5 = e >= e;
  EXPECT_EQ(t5, 1);

  bool t6 = e >= me;
  EXPECT_EQ(t6, 1);

  bool t7 = me >= e;
  EXPECT_EQ(t7, 0);

  // n > m
  BigInt pipi("-314159654");

  bool t8 = pipi >= mpi;
  EXPECT_EQ(t8, 0);

  bool t9 = pi >= e;
  EXPECT_EQ(t9, 1);

  bool t10 = e >= pi;
  EXPECT_EQ(t10, 0);
}

TEST_F(BigIntTest, LT) {
  // n == m
  BigInt e("271828");
  BigInt pi("314159");

  BigInt me("-271828");
  BigInt mpi("-314159");

  bool t1 = e < pi;
  EXPECT_EQ(t1, 1);

  bool t2 = pi < e;
  EXPECT_EQ(t2, 0);

  bool t3 = me < mpi;
  EXPECT_EQ(t3, 0);

  bool t4 = mpi < me;
  EXPECT_EQ(t4, 1);

  bool t5 = e < e;
  EXPECT_EQ(t5, 0);

  bool t6 = e < me;
  EXPECT_EQ(t6, 0);

  bool t7 = me < e;
  EXPECT_EQ(t7, 1);

  // n > m
  BigInt pipi("-314159654");

  bool t8 = pipi < mpi;
  EXPECT_EQ(t8, 1);

  bool t9 = pi < e;
  EXPECT_EQ(t9, 0);

  bool t10 = e < pi;
  EXPECT_EQ(t10, 1);
}

TEST_F(BigIntTest, GT) {
  // n == m
  BigInt e("271828");
  BigInt pi("314159");

  BigInt me("-271828");
  BigInt mpi("-314159");

  bool t1 = e > pi;
  EXPECT_EQ(t1, 0);

  bool t2 = pi > e;
  EXPECT_EQ(t2, 1);

  bool t3 = me > mpi;
  EXPECT_EQ(t3, 1);

  bool t4 = mpi > me;
  EXPECT_EQ(t4, 0);

  bool t5 = e > e;
  EXPECT_EQ(t5, 0);

  bool t6 = e > me;
  EXPECT_EQ(t6, 1);

  bool t7 = me > e;
  EXPECT_EQ(t7, 0);

  // n > m
  BigInt pipi("-314159654");

  bool t8 = pipi > mpi;
  EXPECT_EQ(t8, 0);

  bool t9 = pi > e;
  EXPECT_EQ(t9, 1);

  bool t10 = e > pi;
  EXPECT_EQ(t10, 0);
}

TEST_F(BigIntTest, EqualityTests) {
  EXPECT_EQ(e, e);
  EXPECT_NE(e, pi);
}

TEST_F(BigIntTest, NotEqual) {

  BigInt pi("314159");
  BigInt e("271828");

  bool t1 = e != pi;
  EXPECT_EQ(t1, 1);

  bool t2 = e != e;
  EXPECT_EQ(t2, 0);
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
  EXPECT_EQ(c2.xleft, nan);
  EXPECT_EQ(c2.xright, nan);
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
  BigInt eight("8");
  BigInt three("3");
  BigInt mthree("-3");

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

  BigInt truth8("5");
  EXPECT_EQ(eight - three, truth8);

  BigInt truth9("11");
  EXPECT_EQ(eight - mthree, truth9);
}

TEST_F(BigIntTest, PreIncrementTest) {

  BigInt z("234324324234");
  ++z;
  BigInt truth("234324324235");
  EXPECT_EQ(z, truth);
}

TEST_F(BigIntTest, PreDecrementTest) {
  BigInt z("234324324234");
  --z;
  BigInt truth("234324324233");
  EXPECT_EQ(z, truth);
}

TEST_F(BigIntTest, GetNumerusPtrTest) {

  BigInt z("234324324234");
  std::unique_ptr<std::vector<uint8_t>> v = z.numerus_ptr();
  std::vector<uint8_t> v1 = {2, 3, 4, 3, 2, 4, 3, 2, 4, 2, 3, 4};
  EXPECT_EQ(*v, v1);
}

int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
