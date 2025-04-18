import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.Random;

public class testgenerator2
{
    public static void main(String args[]) throws IOException
    {
        if(args.length < 4)
        {
            System.out.println("Usage java testgenerator [INPUT] [OUTPUT] [NUM TEST CASE] [MAX NUM BITS]");
            return;
        }
        BufferedReader br = new BufferedReader (new InputStreamReader(System.in));

        FileWriter out = new FileWriter(args[0]);
        FileWriter ans = new FileWriter(args[1]);

        int DATA = Integer.parseInt(args[2]);
        int maxNumBits = Integer.parseInt(args[3]);

        String[] operators = {"+", "-", "*", "/", "%"};
        Random rnd = new Random();

        System.out.println("Generating " + DATA + " data.");
        System.out.println("Input " + args[0]);
        System.out.println("Output " + args[1]);
        System.out.println("Number of test cases " + DATA);
        System.out.println("Max number of bits " + maxNumBits);

        for(int i = 0; i < DATA; i++)
        {
            System.out.println( (i+1) + "...");
            BigInteger a = new BigInteger (maxNumBits, rnd);
            BigInteger b = new BigInteger (maxNumBits - rnd.nextInt(maxNumBits) , rnd);

            String op = operators[rnd.nextInt(operators.length)];
            System.out.println("Operator " + op);
            
            out.write(a.toString() + ' ' + op + ' ' + b.toString());

            BigInteger c = new BigInteger("0");
            switch(op)
            {
                case "+":
                    c = a.add(b);
                    break;
                case "-":
                    c = a.subtract(b);
                    break;
                case "*":
                    c = a.multiply(b);
                    break;
                case "/":
                    c = a.divide(b);
                    break;
                case "%":
                    c = a.mod(b);
                    break;
            }
            ans.write(c.toString());

            if(i < DATA - 1)
            {
                ans.write('\n');
                out.write('\n');
            }
            out.flush();
            ans.flush();
        }
        out.close();
        ans.close();
        System.out.println("Done.");
    }
}
