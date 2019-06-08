import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.Random;

public class PerfTestGenerator
{
    public static void main(String args[]) throws IOException
    {
        if(args.length < 4)
        {
            System.out.println("Usage java perftestgenerator [INPUT] [OUTPUT] [NUM DATA] [BIT]");
            return;
        }
        BufferedReader br = new BufferedReader (new InputStreamReader(System.in));

        FileWriter out = new FileWriter(args[0]);
        FileWriter ans = new FileWriter(args[1]);
        int DATA = Integer.parseInt(args[2]);

        Random rnd = new Random();

        int bitIncrement = Integer.parseInt(args[3]);

        int maxNumBits = bitIncrement;

        for(int i = 0; i < DATA; i++)
        {
            System.out.println( (i+1) + ": " + maxNumBits);
            BigInteger a = new BigInteger (maxNumBits, rnd);
            BigInteger b = new BigInteger (maxNumBits, rnd);
            
            out.write(a.toString() + " * " + b.toString());
            // ans.write("Data : " + (i+1) + ' ');
            BigInteger c = a.multiply(b);;
            ans.write(c.toString());

            if(i < DATA - 1)
            {
                ans.write('\n');
                out.write('\n');
            }
            out.flush();
            ans.flush();

            maxNumBits += bitIncrement;
        }
        out.close();
        ans.close();
        System.out.println("Done.");
    }
}
