package com.hartwig.hmftools.neo.bind;

import static java.lang.Math.pow;

import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_AMINO_ACID;
import static com.hartwig.hmftools.neo.bind.BindConstants.aminoAcidIndex;
import static com.hartwig.hmftools.neo.bind.BindData.DELIM;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;

public class BlosumMapping
{
    // mappings from one letter to another, indexed by the standard order so consistent with other data structures which reference AAs
    private final int[][] mMappings;

    public BlosumMapping()
    {
        int aminoAcidCount = BindConstants.AMINO_ACIDS.size();
        mMappings = new int[aminoAcidCount][aminoAcidCount];
        load();
    }

    public int selfMapping(final int aaIndex)
    {
        return mMappings[aaIndex][aaIndex];
    }

    public int selfMapping(final char aa)
    {
        int aaIndex = aminoAcidIndex(aa);
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
        //if(aa1Index == INVALID_AMINO_ACID || aa2Index == INVALID_AMINO_ACID)
        //    return 0;

        return mMappings[aa1Index][aa2Index];
    }

    public double calcSequenceBlosumScore(final String sequence)
    {
        int total = 0;

        for(int i = 0; i < sequence.length(); ++i)
        {
            char aminoAcid = sequence.charAt(i);
            int mapping = selfMapping(aminoAcid);
            total += mapping;
        }

        return pow(2, total);
    }

    public double calcSequenceBlosumScore(final String seq1, final String seq2)
    {
        if(seq1.length() != seq2.length())
            return 0;

        if(seq1.equals(seq2))
            return calcSequenceBlosumScore(seq1);

        double total = 0;

        for(int i = 0; i < seq1.length(); ++i)
        {
            char aa1 = seq1.charAt(i);
            char aa2 = seq2.charAt(i);
            int mapping = map(aa1, aa2);
            total += mapping;
        }

        return pow(2, total);
    }

    private void load()
    {
        final List<String> lines = new BufferedReader(new InputStreamReader(
                RefGenomeCoordinates.class.getResourceAsStream("/ref/blosum62.csv")))
                .lines().collect(Collectors.toList());

        String[] columns = lines.get(0).split(DELIM);
        lines.remove(0);

        if(columns.length != 21)
            return;

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
            String[] items = line.split(DELIM);

            if(items.length != 21)
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
