package com.hartwig.hmftools.linx.utils;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.neoepitope.RefGenomeInterface;

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

        if(chrBases != null && chrBases.length() > posEnd)
            return chrBases.substring((int)posStart, (int)posEnd + 1);
        else
            return "";
    }

}
