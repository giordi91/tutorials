#include <gmock/gmock.h>

class myClass
{
public:
	myClass(int value):m_value(value){};
	int m_value;
};


TEST(myClass, constructor_test)
{
	myClass c(10);
	EXPECT_EQ(c.m_value, 10 );
}


TEST(myClass, constructor_test2)
{
	myClass c(15);
	EXPECT_EQ(c.m_value, 15 );
}