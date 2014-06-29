import java.io.*;
import java.math.*;
import java.util.*;

public class chk_primes {
  public static void main(String[] args) throws Exception {
    long start = System.currentTimeMillis();
    Scanner s = new Scanner(System.in);
    for (;;) {
      BigInteger n;
      try {
        n = new BigInteger(s.next());
      } catch (NoSuchElementException e) {
        break;
      }
      if (n.isProbablePrime(10)) {
        System.out.println("YES");
      } else {
        System.out.println("NO");
      }
    }
    long end = System.currentTimeMillis();
    System.out.println(end - start);
  }
}
