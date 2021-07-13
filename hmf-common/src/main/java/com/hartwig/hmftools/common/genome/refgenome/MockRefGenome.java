package com.hartwig.hmftools.common.genome.refgenome;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

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
                .filter(x -> x[0] >= 0 && x[1] < chrBases.length())
                .forEach(x -> chrBases.substring(x[0], x[1] + 1));

        return refBases.toString();
    }

}
