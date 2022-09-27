package com.hartwig.hmftools.cider;

import static java.lang.Math.pow;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.base.Preconditions;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.codon.Codons;

public class BlosumMapping
{
    public static final byte INVALID_AMINO_ACID = -1;
    private static final int DEFAULT_STOP_CODON_PENALTY = -10;

    public static final char[] AMINO_ACIDS = {
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'};

    // an array that allows fast conversion from the ascii code of the amino acid to index in the above AMINO_ACIDS list
    private static final byte[] AMINO_ACID_INDEX = new byte[charToInt('Y') + 1];

    // mappings from one letter to another.
    private final byte[] mMappings = new byte[calcMappingArrayLength()];

    private static final String BLOSUM62_MATRIX =
              "AminoAcid,A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V\n"
            + "A,4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0\n"
            + "R,-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3\n"
            + "N,-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3\n"
            + "D,-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3\n"
            + "C,0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1\n"
            + "Q,-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2\n"
            + "E,-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2\n"
            + "G,0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3\n"
            + "H,-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3\n"
            + "I,-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3\n"
            + "L,-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1\n"
            + "K,-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2\n"
            + "M,-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1\n"
            + "F,-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1\n"
            + "P,-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2\n"
            + "S,1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2\n"
            + "T,0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0\n"
            + "W,-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3\n"
            + "Y,-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1\n"
            + "V,0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4";

    static
    {
        Arrays.fill(AMINO_ACID_INDEX, INVALID_AMINO_ACID);

        // assign an index for each amino acid
        for(byte i = 0; i < AMINO_ACIDS.length; ++i)
        {
            AMINO_ACID_INDEX[charToInt(AMINO_ACIDS[i])] = i;
        }
    }

    private static int charToInt(char c) { return Character.toUpperCase(c) - 'A'; }

    public static int aminoAcidIndex(final char aminoAcid)
    {
        int charInt = charToInt(aminoAcid);
        if (charInt < 0 || charInt >= AMINO_ACID_INDEX.length)
            return -1;
        return AMINO_ACID_INDEX[charInt];
    }

    public BlosumMapping(int stopCodonPenalty)
    {
        Preconditions.checkArgument(stopCodonPenalty >= Byte.MIN_VALUE);
        Preconditions.checkArgument(stopCodonPenalty <= 0);
        load((byte)stopCodonPenalty);
    }

    public BlosumMapping()
    {
        this(DEFAULT_STOP_CODON_PENALTY);
    }

    private int getMapping(int aa1Index, int aa2Index)
    {
        if (aa1Index > aa2Index)
        {
            int tmp = aa1Index;
            aa1Index = aa2Index;
            aa2Index = tmp;
        }
        return mMappings[getMappingArrayIndex(aa1Index, aa2Index)];
    }

    private void setMapping(int aa1Index, int aa2Index, byte val)
    {
        mMappings[getMappingArrayIndex(aa1Index, aa2Index)] = val;
    }

    // we store only half of the matrix, so each row is smaller than previous
    // i.e. for the n x n where n = 4 matrix:
    //             row index(i)    row array index        check (i + 1) * i / 2
    //    0            0              0 = 0                 (0 + 1) * 0 / 2 = 0
    //    1 2          1              1 = 0 + 1             (1 + 1) * 1 / 2 = 1
    //    3 4 5        2              3 = 0 + 1 + 2         (2 + 1) * 2 / 2 = 3
    //    6 7 8 9      3              6 = 0 + 1 + 2 + 3     (3 + 1) * 3 / 2 = 6
    //
    // from here we can work out that for row index i,
    // row array index is given by arithmetic sum:
    // S[0 .. i] = (i + 1) * i / 2
    // we can see in the check column that the values do match
    private static int getMappingArrayIndex(int row, int col)
    {
        if (col > row)
        {
            // swap them around
            int tmp = row;
            row = col;
            col = tmp;
        }

        Preconditions.checkArgument(row >= col);
        int rowArrayIndex = (row + 1) * row / 2;
        return rowArrayIndex + col;
    }

    private int calcMappingArrayLength()
    {
        // using N as row index we get the index of the last element + 1
        // which is the array length
        return getMappingArrayIndex(AMINO_ACIDS.length, 0);
    }

    public int selfMapping(final char aa)
    {
        int aaIndex = aminoAcidIndex(aa);

        if (aaIndex == INVALID_AMINO_ACID)
            throw new IllegalArgumentException("invalid amino acid: " + aa);

        return getMapping(aaIndex, aaIndex);
    }

    public int map(final char aa1, final char aa2)
    {
        int aa1Index = aminoAcidIndex(aa1);
        int aa2Index = aminoAcidIndex(aa2);

        // better to fail than return a potentially incorrect value
        if (aa1Index == INVALID_AMINO_ACID)
            throw new IllegalArgumentException("invalid amino acid: " + aa1);

        if (aa2Index == INVALID_AMINO_ACID)
            throw new IllegalArgumentException("invalid amino acid: " + aa2);

        return getMapping(aa1Index, aa2Index);
    }

    public int calcSequenceSum(final String sequence)
    {
        int total = 0;

        for(int i = 0; i < sequence.length(); ++i)
        {
            char aminoAcid = sequence.charAt(i);
            int mapping = selfMapping(aminoAcid);
            total += mapping;
        }

        return total;
    }

    public int calcSequenceSum(final String seq1, final String seq2)
    {
        if(seq1.length() != seq2.length())
            return 0;

        if(seq1.equals(seq2))
            return calcSequenceSum(seq1);

        int total = 0;

        for(int i = 0; i < seq1.length(); ++i)
        {
            char aa1 = seq1.charAt(i);
            char aa2 = seq2.charAt(i);
            int mapping = map(aa1, aa2);
            total += mapping;
        }

        return total;
    }

    public double calcSequenceScore(final String sequence)
    {
        return pow(2, calcSequenceSum(sequence));
    }

    public double calcSequenceScore(final String seq1, final String seq2)
    {
        return pow(2, calcSequenceSum(seq1, seq2));
    }

    public double calcSequenceSimilarity(final String seq1, final String seq2)
    {
        if(seq1.length() != seq2.length())
            return -1;

        if(seq1.equals(seq2))
            return 0;

        double similarity = 0;

        for(int i = 0; i < seq1.length(); ++i)
        {
            char aa1 = seq1.charAt(i);
            char aa2 = seq2.charAt(i);

            int bs1 = selfMapping(aa1);
            int bs2 = selfMapping(aa2);
            int map = map(aa1, aa2);

            similarity += (bs1 + bs2) * 0.5 - map;
        }

        return similarity;
    }

    private void load(byte stopCodonPenalty)
    {
        final List<String> lines = new ArrayList<>(Arrays.asList(BLOSUM62_MATRIX.split("\n")));

        String[] columns = lines.get(0).split(",");
        lines.remove(0);

        if(columns.length != AMINO_ACIDS.length)
        {
            throw new RuntimeException("invalid blosum62 input file, number of column != 21");
        }

        Map<Integer,Integer> columnAaMap = Maps.newHashMap();

        for(int i = 1; i < columns.length; ++i)
        {
            char aa = columns[i].charAt(0);
            int aaIndex = aminoAcidIndex(aa);

            if(aaIndex == INVALID_AMINO_ACID)
                return;

            columnAaMap.put(i, aaIndex);
        }

        for(String line : lines)
        {
            String[] items = line.split(",");

            if(items.length != AMINO_ACIDS.length)
                return;

            char aa1 = items[0].charAt(0);
            int aa1index = aminoAcidIndex(aa1);

            if(aa1index == INVALID_AMINO_ACID)
                return;

            for(int i = 1; i < items.length; ++i)
            {
                byte correlation = Byte.parseByte(items[i]);
                int aa2index = columnAaMap.get(i);

                if(aa2index == INVALID_AMINO_ACID)
                    return;

                setMapping(aa1index, aa2index, correlation);
            }
        }

        // fill in the stop codon rows
        int stopIndex = aminoAcidIndex(Codons.STOP_AMINO_ACID);
        for (char aa : AMINO_ACIDS)
        {
            int aaIndex = aminoAcidIndex(aa);
            setMapping(aaIndex, stopIndex, stopCodonPenalty);
            setMapping(stopIndex, aaIndex, stopCodonPenalty);
        }

        // no penalty for stop codon self mapping
        setMapping(stopIndex, stopIndex, (byte)0);
    }
}
