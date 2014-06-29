import java.io.*;
import java.math.*;
import java.util.*;

public class gen_primes {
  public static void main(String[] args) throws Exception {
    Random r = new Random();

    final int TESTN = 100;
    final int BASE = (int)4e17;
    for (int i = 0; i < TESTN; i++) {
      BigInteger a = new BigInteger(59, r).add(BigInteger.valueOf(BASE)).nextProbablePrime();
      System.out.println(a);
    }
  }
}
