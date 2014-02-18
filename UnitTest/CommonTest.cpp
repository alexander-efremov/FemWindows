#include "gtest/gtest.h"

inline double round (double x)
{
	return ceil (x + 0.5);
}


TEST(CommonTest, SplitIntoChunksTest)
{
	int count = 9;

	int data[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};

	int chunk = 5;

	int parts = round(9/chunk);
	for (int i = 0; i < parts; ++i)
	{
		int printf1 = printf("Chunk ");
		int printf2 = printf("%d\n", i);
		int offset = i*chunk;
		for (int j = 0; j < chunk && (j + offset) < count; ++j)
		{
			printf("%d\t", data[j + offset]);
		}
		printf("\n");
	}
}