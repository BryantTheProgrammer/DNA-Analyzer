#include "pch.h"
//#include "gtest/gtest.h"
#include "DNAStrand.h"

TEST(DNAStrand, TestConstants) {
  DNAStrand cS1 = { ">S1","AAGC" };
  EXPECT_EQ(cS1.getBase(1), cS1.getBase(0));
  ASSERT_EQ ('C', cS1.getBase (3));
  EXPECT_EQ ('G', cS1.getBase (2));
  EXPECT_TRUE(true);
}