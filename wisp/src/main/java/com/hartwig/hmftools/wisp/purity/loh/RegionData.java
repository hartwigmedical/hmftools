package com.hartwig.hmftools.wisp.purity.loh;

import static java.lang.Math.max;

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
}
