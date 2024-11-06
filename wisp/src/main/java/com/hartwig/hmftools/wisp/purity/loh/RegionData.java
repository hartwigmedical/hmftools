package com.hartwig.hmftools.wisp.purity.loh;

import static java.lang.Math.max;
import static java.lang.Math.round;
import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;

public class RegionData
{
    public final PurpleCopyNumber CopyNumber;
    public final List<SiteData> Sites;

    public RegionData(final PurpleCopyNumber copyNumber)
    {
        CopyNumber = copyNumber;
        Sites = Lists.newArrayList();
    }

    public double averageAF()
    {
        double totalAF = Sites.stream().mapToDouble(x -> x.AF).sum();
        return Sites.size() > 0 ? totalAF / Sites.size() : 0;
    }

    public double impliedPurity()
    {
        double avgAf = averageAF();
        double denom = avgAf * (CopyNumber.averageTumorCopyNumber() - 2) + 1;
        return denom > 0 ? max((1 - 2 * avgAf) / denom, 0) : 0;
    }

    public int expectedSupport(double purity)
    {
        double denom = purity * (CopyNumber.averageTumorCopyNumber() - 2) + 2;

        if(denom <= 0)
            return 0;

        double expectedAf = (1 - purity) / denom;
        int expectedSupport = (int)round(expectedAf * siteDepth());
        return expectedSupport;
    }

    public int siteSupport() { return Sites.stream().mapToInt(x -> x.Support).sum(); }
    public int siteDepth() { return Sites.stream().mapToInt(x -> x.SampleDepth).sum(); }

    public String toString()
    {
        return format("region(%s:%d-%d) site(%d) af(%.3f) purity(%.5f) frags(depth=%d support=%d)",
                CopyNumber.chromosome(), CopyNumber.start(), CopyNumber.end(), Sites.size(),
                averageAF(), impliedPurity(), siteDepth(), siteSupport());
    }
}
