import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.Random;

public class DivPerfTestGenerator
{
    public static void main(String args[]) throws IOException
    {
        if(args.length < 4)
        {
            System.out.println("Usage java DivPerfTestGenerator [INPUT] [OUTPUT] [NUM DATA] [BIT]");
            return;
        }
        BufferedReader br = new BufferedReader (new InputStreamReader(System.in));

        FileWriter out = new FileWriter(args[0]);
        FileWriter ans = new FileWriter(args[1]);
        int DATA = Integer.parseInt(args[2]);

        Random rnd = new Random();

        int bitIncrement = Integer.parseInt(args[3]);

        // we want the divisor to be smaller than the dividend
        int maxNumBitsA = bitIncrement;
        int maxNumBitsB = maxNumBitsA/2;

        for(int i = 0; i < DATA; i++)
        {
            System.out.println( (i+1) + ": " + maxNumBitsA + " " + maxNumBitsB);
            BigInteger a = new BigInteger (maxNumBitsA, rnd);
            BigInteger b = new BigInteger (maxNumBitsB, rnd);
            
            out.write(a.toString() + " / " + b.toString());
            // ans.write("Data : " + (i+1) + ' ');
            BigInteger[] c = a.divideAndRemainder(b);;
            ans.write(c[0].toString());
            ans.write('\n');
            ans.write(c[1].toString());

            if(i < DATA - 1)
            {
                ans.write('\n');
                out.write('\n');
            }
            out.flush();
            ans.flush();

            maxNumBitsA += bitIncrement;
            maxNumBitsB = maxNumBitsA/2;
        }
        out.close();
        ans.close();
        System.out.println("Done.");
    }
}
