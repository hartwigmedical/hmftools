package com.hartwig.hmftools.common.purple.region;

import java.util.List;

import com.hartwig.hmftools.common.purple.FittedCopyNumber;

public class ConsolidatedRegionFactory {
    public static List<ConsolidatedRegion> broad(List<FittedCopyNumber> copyNumbers) {
        return new BroadRegions().broad(copyNumbers);
    }

    public static List<ConsolidatedRegion> smooth(List<FittedCopyNumber> copyNumbers, List<ConsolidatedRegion> broadRegions) {
        return new SmoothedRegions(copyNumbers).smooth(broadRegions);
    }

}
