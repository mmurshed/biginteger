#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>

using namespace std;

#include "../BigInteger.h"
#include "../BigIntegerIO.h"
#include "ExpressionEvaluator.h"

using namespace BigMath;

int main()
{
	cout << "Enter an expression (enter to exit)" << endl;
	while (true)
	{
		// Get a string from console
		string line;
		std::getline(cin, line);

		// If the string is empty, then exit
		if (line.size() == 0)
			break;

		// Evaluate the expression
		ExpressionEvaluator eval;
		try
		{
			BigInteger res = eval.Eval(line.c_str());
			cout << res << endl;
		}
		catch (const char *msg)
		{
			cerr << msg << endl;
		}
	}
	return 0;
}
