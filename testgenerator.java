import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.Random;

public class testgenerator
{

    public static void main(String args[])
    throws IOException
    {
        BufferedReader br = new BufferedReader (new InputStreamReader(System.in));

        System.out.print("Enter number of DATA to generate : ");
        int DATA = Integer.parseInt(br.readLine());

        FileWriter out = new FileWriter("test.in");
        FileWriter ans = new FileWriter("test.ans");

        System.out.println("Generating " + DATA + " data.\nPlease wait...");

        Random rnd = new Random();

        for(int i=0;i<DATA;i++)
        {
            int numBits = 16600;
            System.out.println( (i+1) + "...");
            BigInteger a = new BigInteger (numBits, rnd);
            BigInteger b = new BigInteger (numBits - rnd.nextInt(numBits) , rnd);
            out.write(a.toString() + " / " + b.toString());
            ans.write("Data : " + (i+1) + " " + a.divide(b).toString() + "\n");
            out.write('\n');
            out.flush();
        }
        out.close();
        ans.close();
        System.out.println("Done.");
    }
}
