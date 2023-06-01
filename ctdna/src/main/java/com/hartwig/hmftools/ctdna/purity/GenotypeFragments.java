package com.hartwig.hmftools.ctdna.purity;

import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNT;
import static com.hartwig.hmftools.common.variant.SageVcfTags.UMI_TYPE_COUNTS;

public class GenotypeFragments
{
    public final String SampleName;
    public final int AlleleCount;
    public final int Depth;
    public final double QualTotal;
    public final int[] UmiTypeCounts;

    public GenotypeFragments(final String sampleName, final int alleleCount, final int depth, final double qualTotal, final Object umiTypeCountsRaw)
    {
        SampleName = sampleName;
        AlleleCount = alleleCount;
        Depth = depth;
        QualTotal = qualTotal;

        UmiTypeCounts = new int[UMI_TYPE_COUNT];

        if(umiTypeCountsRaw != null)
        {
            String[] umiValues = ((String)umiTypeCountsRaw).split(",", UMI_TYPE_COUNT);

            if(umiValues.length == UMI_TYPE_COUNT)
            {
                for(int i = 0; i < UmiTypeCounts.length; ++i)
                {
                    UmiTypeCounts[i] = Integer.parseInt(umiValues[i]);
                }
            }
        }
    }

    public double qualPerAlleleFragment()
    {
        return AlleleCount > 0 ? QualTotal / AlleleCount : 0;
    }
}
