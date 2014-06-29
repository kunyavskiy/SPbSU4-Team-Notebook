import java.io.*;
import java.math.*;
import java.util.*;

public class gen_mult {
  public static void main(String[] args) throws Exception {
    Random r = new Random();

    final int TESTN = 100;
    final int BASE = (int)4e8;
    for (int i = 0; i < TESTN; i++) {
      BigInteger a = new BigInteger(29, r).add(BigInteger.valueOf(BASE)).nextProbablePrime();
      BigInteger b = new BigInteger(29, r).add(BigInteger.valueOf(BASE)).nextProbablePrime();
      System.out.println(a.multiply(b));
    }
  }
}
