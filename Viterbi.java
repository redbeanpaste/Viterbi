package viterbi;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;

public class Viterbi {

    Object[] readInput() throws FileNotFoundException, IOException
    {
        BufferedReader br = new BufferedReader(new FileReader("states.txt"));
        int state = 0;
        String line; 
        // read state
        while(br.readLine() != null) state += 1;
        double[][] tranTable = new double[state+2][state+2];
        String[] states = new String[state+1]; int i = 1;
        states[0] = "0";
        br = new BufferedReader(new FileReader("states.txt"));
        while((line = br.readLine()) != null) 
        {
            states[i] = line;
            i += 1;
        }
        br = new BufferedReader(new FileReader("symbols.txt"));
        // read symbol
        String sym = "";
        while((line = br.readLine()) != null) 
        {
            sym = sym.concat(line);
        }
        String [] symbol = sym.split("");
        double[][] emisTable = new double[state+2][symbol.length];
        br = new BufferedReader(new FileReader("starting.txt"));
        // read startprob
        int j = 1;
        while((line = br.readLine()) != null) 
        {
            tranTable[0][j] = Double.valueOf(line);
            j += 1;
        }
        br = new BufferedReader(new FileReader("transition.txt"));
        int row = 1;
        double sum = 0, value;
        // read tranprob
        while((line = br.readLine()) != null) 
        {
            String[] read = line.split(" ");
            for(int col = 1; col <= read.length; col++)
            {
                value = Double.valueOf(read[col-1]);
                tranTable[row][col] = value;
                sum += value;
            }
            if( sum < 1.0 ) tranTable[row][tranTable.length-1] = round((1.0-sum),1); 
            sum = 0;
            row += 1;
        }
        br = new BufferedReader(new FileReader("emission.txt"));
        row = 1;
        // read emisprob
        while((line = br.readLine()) != null) 
        {
            String[] read = line.split(" ");
            for(int col = 0; col < read.length; col++)
            {
                emisTable[row][col] = Double.valueOf(read[col]);
            } 
            row += 1;
        }
        br = new BufferedReader(new FileReader("sequence.txt"));
        // read seq
        line = br.readLine();
        String[] seq = line.split("");
        //create path table
        int N = tranTable.length;
        int T = seq.length + 1;
        double[][] pathTable = new double[N][T];
        
        System.out.println("state "+Arrays.toString(states));
        System.out.println("tran");
        printMatrix(tranTable);
        System.out.println("emis");
        printMatrix(emisTable);
        
        Object[] pkg = {tranTable, emisTable, pathTable, symbol, seq, states};
        return pkg;
    }
    
    Object[] calcTable(Object[] obj)
    {
        double[][] tranTable = (double[][]) obj[0];
        double[][] emisTable = (double[][]) obj[1];
        double[][] pathTable = (double[][]) obj[2];
        String[] symbol = (String[]) obj[3];
        String[] seq = (String[]) obj[4];
        int N = tranTable.length;
        int T = seq.length;
        int[][] ptr = new int[N][T+1];
        init(pathTable, seq.length);
        
        System.out.println("seq: "+Arrays.toString(seq));
        
        // calculate begin
        for(int t = 1; t <= T; t++)
        {
            String sym = seq[t-1];
            for(int k = 1; k < N; k++)
            {
                double e = getEmisProb(emisTable, symbol, sym, k);
                double maxProb = 0;
                double prob;
                for(int l = 1; l < N ; l++)
                {
                    prob = e*tranTable[l-1][k]*pathTable[l-1][t-1];
                    if(prob > maxProb)
                    {
                        maxProb = round(prob,20);
                        ptr[k][t] = l;
                    }
                }
                pathTable[k][t] = maxProb;
            }
        }
        
        String[] states = (String[]) obj[5];
        Object[] pkg = {tranTable, emisTable, pathTable, ptr, symbol, seq, states};
        termination(pkg);
        
        System.out.println("path");
        printMatrix(pathTable);
        System.out.println("ptr");
        printMatrix(ptr);
        
        return pkg;
        
    }
    
    void init(double[][] pathTable, int col)
    {
        pathTable[0][0] = 1;
        for(int row = 1; row < pathTable.length; row++) pathTable[row][0] = 0;
        for(int len = 1; len < col; len++) pathTable[0][len] = 0;
    }
    
    double getEmisProb(double[][] emisTable, String[] symbol, String e, int state)
    {
        int sym;
        for(sym = 0; sym < symbol.length; sym++)
        {
            if(symbol[sym].equals(e)) 
                break;
        }
        return emisTable[state][sym];
    }
    
    void termination(Object[] obj)
    {
        double[][] tranTable = (double[][]) obj[0];
        double[][] emisTable = (double[][]) obj[1];
        double[][] pathTable = (double[][]) obj[2];
        int[][] ptr = (int[][]) obj[3];
        String[] seq = (String[]) obj[5];
        
        int t = seq.length;
        int N = tranTable.length;
        int probPath = 0;
        double maxProb = 0, prob;
        for(int k = 0; k < N; k++)
        {
            prob = pathTable[k][t]*tranTable[k][N-1];
            if(prob > maxProb)
            {
                maxProb = round(prob,20);
                probPath = k;
            }
        }
        pathTable[N-1][t] = maxProb;
        ptr[N-1][t] = probPath + 1;
        
    }
    
    void traceback(Object[] obj)
    {
        double[][] pathTable = (double[][]) obj[2];
        int[][] ptr = (int[][]) obj[3];
        String[] seq = (String[]) obj[5];
        String[] states = (String[]) obj[6];
        int T = seq.length;
        int N = ptr.length;
        int k = N-1;
        int[] trb = new int[T+1];
        int probPath;
        trb[T] = ptr[k][T];
        for(int t = T; t > 1; t--)
        {
            probPath = ptr[k-1][t];
            k = probPath;
            trb[t-1] = probPath;
        }
        
        
        System.out.print("Sequence: ");
        for(String s: seq) System.out.print(s);
        System.out.print("\nP(");
        for(String s: seq) System.out.print(s);
        System.out.print(", ");
        for(int i: trb) 
        {
            if(i == 0) System.out.print("0 ");
            else
            {
                if(i == 1) System.out.print(states[i]+" ");
                else System.out.print(states[i-1]+" ");
            }
        }
        System.out.println(") = "+pathTable[N-1][T]);
        
    }
    
    void findPath() throws IOException
    {
        Object[] pkg = readInput();
        pkg = calcTable(pkg);
        traceback(pkg);
    }
    
    public static void main(String[] args) throws IOException {
        // there is 6 input file for this viterbi algorithm
        // states.txt contains all states of model(start and end state dont include here)
        // symbols.txt contains all symbols that present in sequence and emitted from all states
        // starting.txt contains prob from starting state to all states in model
        // transition.txt contains prob from state i to state j in model
        // emission.txt contains prob of state j will emit symbol i
        // sequence.txt contains sequence to find most probably path from model
    
        Viterbi v = new Viterbi();
        v.findPath();
    }
    
    double max(double x, double y, double z)
    {
        if( x >= y && x >= z ) return x;
        else if( y >= x && y >= z ) return y;
        else return z;
    }
    
    void printMatrix(double[][] matrix)
    {
        for(double[] row: matrix)
        {
            for(double num: row)
            {
                System.out.print(num+"\t");
            }
            System.out.println("");
        }
        System.out.println("");
    }
    
    void printMatrix(int[][] matrix)
    {
        for(int[] row: matrix)
        {
            for(int read: row)
            {
                System.out.print(read+"\t");
            }
            System.out.println("");
        }
        System.out.println("");
    }
    
    double round(double value, int places)
    {
        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }
    
    /*
    matrix in needed ; n is no of states
        - transitionTable (n+2)*(n+2)
        - emissionTable n*4(ACGT)
        - pathTable (n+2)*(seq_len+1)
    */
}
