/**
 * Version 9.0
 * S. M. Mahbub Murshed (murshed@gmail.com)
 */

/*
 * From https://www.strchr.com/expression_evaluator
 * 
 * Recursive descent parsing
 * 
 * A popular approach (called recursive descent parsing) views 
 * an expression as a hierarchical structure. On the top level,
 * an expression consists of several summands:
 * 
 *  1.5  +  2  +  3.4 * (25 – 4) / 2  –  8
 * 
 * [1.5] + [2] + [3.4 * (25 – 4) / 2] – [8]
 * 
 * Summands are 1.5, 2, 3.4 * (25 – 4) / 2, and 8
 * 
 * Summands are separated by "+" and "–" signs. They include 
 * single numbers (1.5, 2, 8) and more complex expressions 
 * [3.4 * (25 – 4) / 2].
 * 
 * (Strictly speaking, 8 is not a summand here, but a subtrahend.
 * For the purposes of this article, we will refer to the parts 
 * of expression separated by "–" as summands, though a 
 * mathematician would say it's a misnomer.)
 * 
 * Each summand, in turn, consist of several factors (we will 
 * call "factors" the things separated by "*" and "/"):
 * 
 *  3.4  *  (25 – 4)  /  2
 * 
 * [3.4] * [(25 – 4)] / [2]
 * 
 * Factors are 3.4, (25 – 4), and 2. In '1.5', there is one factor, 1.5
 * 
 * The summand "1.5" can be viewed as a degenerated case: it constists 
 * of one factor, 1.5. The summand "3.4 * (25 – 4) / 2" consists of 
 * three factors.
 * 
 * There are 2 types of factors: atoms and subexpressions in 
 * parenthesis. In our simple expression evaluator, an atom is 
 * just a number (in more complex translators, atoms can also 
 * include variable names). Subexpressions can be parsed in the
 *  same way as the whole expression (you again found the summands
 *  inside the brackets, then factors, and then atoms).
 * 
 * The hierarchy of summands and factors
 * 
 * [1.5] + [2] + [3.4  * ( 25  –  4)  /  2] – [8]
 *   ↓      ↓      ↓           ↓         ↓     ↓
 * [1.5] + [2] + [3.4] * [(25  –  4)] / [2] – [8]
 *                             ↓
 *                       ([25  –  4])
 *                         ↓      ↓
 *                        [25] – [4]
 * So, we have an hierarchy, where atoms constitute the lowest 
 * level, then they are combined into factors, and the factors 
 * are combined into summands. A natural way to parse such 
 * hierarchical expression is to use recursion.
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
	int parenthesis;

	// Parse a number or an expression in parenthesis
	BigInteger ParseAtom(const char*& expr)
	{
		// Skip spaces
		while(*expr == ' ')
			expr++;

		// Sign before parenthesis or number
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

		// Parenthesis
		if(*expr == '(')
		{
			expr++;
			parenthesis++;

			BigInteger res = ParseSummands(expr);
			res.SetSign(negative);

			if(*expr != ')')
			{
				throw "Unmatched opening parenthesis";
			}

			expr++;
			parenthesis--;

			return res;
		}

		// It should be a number
		if(*expr < '0' || *expr > '9')
		{
			throw "Invalid character";
		}

		int count = 0;
		BigInteger res = BigIntegerParser::Parse(expr, &count);
		res.SetSign(negative);

		// Advance the pointer and return the result
		expr += count;

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
		parenthesis = 0;
		BigInteger res = ParseSummands(expr);
		// Expression should end, and parenthesis should be zero

		if(parenthesis != 0 || *expr == ')')
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