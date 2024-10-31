package com.hartwig.hmftools.common.genome.gc;

public final class GcCalcs
{
    public static double calcGcPercent(final String sequence)
    {
        int gcCount = 0;
        int atCount = 0;
        for(int i = 0; i < sequence.length(); ++i)
        {
            char base = sequence.charAt(i);
            if(base == 'G' || base == 'C')
                ++gcCount;
            else if(base == 'A' || base == 'T')
                ++atCount;
        }
        int total = gcCount + atCount;
        if(total == 0)
            return 0;
        return gcCount / (double)total;
    }
}
