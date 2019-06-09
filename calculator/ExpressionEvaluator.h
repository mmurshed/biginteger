/**
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

#ifndef EXPRESSION_EVALUATOR
#define EXPRESSION_EVALUATOR

#include "../BigInteger.h"
#include "../BigIntegerParser.h"
#include "../BigIntegerOperations.h"

using namespace BigMath;

class ExpressionEvaluator
{
	private:
	int parenthesis_count;

	// Parse a number or an expression in parenthesis
	BigInteger ParseAtom(const char*& expr)
	{
		// Skip spaces
		while(*expr == ' ')
			expr++;

		// Handle the sign before parenthesis (or before number)
		bool negative = false;

		if(*expr == '-')
		{
			negative = true;
			expr++;
		}

		if(*expr == '+')
		{
			expr++;
		}

		// Check if there is parenthesis
		if(*expr == '(')
		{
			expr++;
			parenthesis_count++;

			BigInteger res = ParseSummands(expr);
			res.SetSign(negative);

			if(*expr != ')')
			{
				throw "Unmatched opening parenthesis";
			}

			expr++;
			parenthesis_count--;

			return res;
		}

		// It should be a number
		if(*expr < '0' || *expr > '9')
		{
			throw "Invalid character";
		}

		int count = 0;
		cout << expr << endl;
		BigInteger res = BigIntegerParser::Parse(expr, &count);
		cout << count << endl;
		cout << res << endl;
		res.SetSign(negative);

		// Advance the pointer and return the result
		expr += count;
		cout << expr << endl;

		return res;
	}

	// Parse multiplication and division
	BigInteger ParseFactors(const char*& expr) {

		BigInteger a = ParseAtom(expr);

		while(true)
		{
			// Skip spaces
			while(*expr == ' ')
				expr++;

			// Save the operation and position
			char op = *expr;
			const char* pos = expr;

			if(op != '/' && op != '*')
				return a;

			expr++;

			BigInteger b = ParseAtom(expr);

			// Perform the saved operation
			if(op == '/')
			{
				// Division by zero
				if(b.IsZero())
				{
					throw "Division by zero";
				}
				a = a / b;
			}
			else
				a = a * b;
		}

		return a;
	}

	// Parse addition and subtraction
	BigInteger ParseSummands(const char*& expr) {
		BigInteger a = ParseFactors(expr);

		while(true)
		{
			// Skip spaces
			while(*expr == ' ')
				expr++;

			char op = *expr;

			if(op != '-' && op != '+')
				return a;

			cout << op << endl;

			expr++;

			BigInteger b = ParseFactors(expr);

			if(op == '-')
				a = a - b;
			else
				a = a + b;
		}

		return a;
	}

public:
	BigInteger Eval(const char* expr) {
		parenthesis_count = 0;
		BigInteger res = ParseSummands(expr);
		// Expression should end, and parenthesis_count should be zero

		if(parenthesis_count != 0 || *expr == ')')
		{
			throw "Parenthesis mismatch";
		}

		if(*expr != '\0')
		{
			throw "Invalid expression";
		}
		return res;
	};
};

#endif