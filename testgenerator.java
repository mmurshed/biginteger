import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.Random;

public class testgenerator
{
    public static void main(String args[]) throws IOException
    {
        if(args.length < 2)
            return;
        BufferedReader br = new BufferedReader (new InputStreamReader(System.in));

        System.out.print("Enter number of DATA to generate: ");
        int DATA = Integer.parseInt(br.readLine());

        System.out.print("Enter max number of bits to generate: ");
        int maxNumBits = Integer.parseInt(br.readLine());

        System.out.print("Enter operator to test: ");
        String op = br.readLine();

        FileWriter out = new FileWriter(args[0]);
        FileWriter ans = new FileWriter(args[1]);

        System.out.println("Generating " + DATA + " data.\nPlease wait...");

        Random rnd = new Random();        

        for(int i = 0; i < DATA; i++)
        {
            System.out.println( (i+1) + "...");
            BigInteger a = new BigInteger (maxNumBits, rnd);
            BigInteger b = new BigInteger (maxNumBits - rnd.nextInt(maxNumBits) , rnd);
            
            out.write(a.toString() + ' ' + op + ' ' + b.toString());
            ans.write("Data : " + (i+1) + ' ');
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
