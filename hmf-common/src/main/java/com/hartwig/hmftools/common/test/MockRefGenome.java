package com.hartwig.hmftools.common.test;

import static com.hartwig.hmftools.common.codon.Nucleotides.DNA_BASES;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

// test implementation of the ref genome
public class MockRefGenome implements RefGenomeInterface
{
    public final Map<String,String> RefGenomeMap;

    public MockRefGenome()
    {
        RefGenomeMap = Maps.newHashMap();
    }

    @Override
    public String getBaseString(final String chromosome, int posStart, int posEnd)
    {
        String chrBases = RefGenomeMap.get(chromosome);

        if(chrBases != null && posStart >= 0 && chrBases.length() > posEnd)
            return chrBases.substring(posStart, posEnd + 1);
        else
            return "";
    }

    @Override
    public String getBaseString(final String chromosome, final List<int[]> baseRanges)
    {
        String chrBases = RefGenomeMap.get(chromosome);

        if(chrBases == null)
            return "";

        StringBuilder refBases = new StringBuilder();

        baseRanges.stream()
                .filter(x -> x[SE_START] >= 0 && x[SE_END] < chrBases.length())
                .forEach(x -> refBases.append(chrBases.substring(x[SE_START], x[SE_END] + 1)));

        return refBases.toString();
    }

    public static String generateRandomBases(int length)
    {
        char[] str = new char[length];

        int baseIndex = 0;
        for(int i = 0; i < length; ++i)
        {
            str[i] = DNA_BASES[baseIndex];

            if(baseIndex == DNA_BASES.length - 1)
                baseIndex = 0;
            else
                ++baseIndex;
        }

        return String.valueOf(str);
    }

    public static String getNextBase(final String base) { return String.valueOf(getNextBase(base.charAt(0))); }

    public static char getNextBase(final char base)
    {
        for(int i = 0; i < DNA_BASES.length; ++i)
        {
            if(DNA_BASES[i] == base)
            {
                return i == DNA_BASES.length - 1 ? DNA_BASES[0] : DNA_BASES[i + 1];
            }
        }

        return base;
    }

}
