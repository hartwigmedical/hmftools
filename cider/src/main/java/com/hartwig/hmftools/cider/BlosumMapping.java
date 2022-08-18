package com.hartwig.hmftools.cider;

import static java.lang.Math.pow;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;

public class BlosumMapping
{
    public static final int INVALID_AMINO_ACID = -1;

    public static final char[] AMINO_ACIDS = {
            'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X'};

    private static final int[] AMINO_ACID_INDEX = new int[charToInt('Y') + 1];

    // mappings from one letter to another, indexed by the standard order so consistent with other data structures which reference AAs
    private final int[][] mMappings;

    static
    {
        Arrays.fill(AMINO_ACID_INDEX, -1);

        // assign an index for each amino acid
        for(int i = 0; i < AMINO_ACIDS.length; ++i)
        {
            AMINO_ACID_INDEX[charToInt(AMINO_ACIDS[i])] = i;
        }
    }

    public static int charToInt(char c) {   return Character.getNumericValue(Character.toUpperCase(c)) - Character.getNumericValue('A'); }

    public static int aminoAcidIndex(final char aminoAcid)
    {
        int charInt = charToInt(aminoAcid);
        if (charInt < 0 || charInt >= AMINO_ACID_INDEX.length)
            return -1;
        return AMINO_ACID_INDEX[charInt];
    }

    public BlosumMapping()
    {
        int aminoAcidCount = AMINO_ACIDS.length;
        mMappings = new int[aminoAcidCount][aminoAcidCount];
        load();
    }

    public int selfMapping(final char aa)
    {
        int aaIndex = aminoAcidIndex(aa);

        if (aaIndex == INVALID_AMINO_ACID)
            throw new IllegalArgumentException("invalid amino acid: " + aa);

        return mMappings[aaIndex][aaIndex];
    }

    public int map(int aa1Index, int aa2Index)
    {
        return mMappings[aa1Index][aa2Index];
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

        return mMappings[aa1Index][aa2Index];
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

    private void load()
    {
        final List<String> lines = new BufferedReader(new InputStreamReader(
                BlosumMapping.class.getClassLoader().getResourceAsStream("blosum62.csv")))
                .lines().collect(Collectors.toList());

        String[] columns = lines.get(0).split(",");
        lines.remove(0);

        if(columns.length != AMINO_ACIDS.length + 1)
        {
            throw new RuntimeException("invalid blosum62 input file, number of column != 22");
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

            if(items.length != AMINO_ACIDS.length + 1)
                return;

            char aa1 = items[0].charAt(0);
            int aa1index = aminoAcidIndex(aa1);

            if(aa1index == INVALID_AMINO_ACID)
                return;

            for(int i = 1; i < items.length; ++i)
            {
                int correlation = Integer.parseInt(items[i]);
                int aa2index = columnAaMap.get(i);

                if(aa2index == INVALID_AMINO_ACID)
                    return;

                mMappings[aa1index][aa2index] = correlation;
            }
        }
    }
}
