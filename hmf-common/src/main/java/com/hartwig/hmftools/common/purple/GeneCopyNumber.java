package com.hartwig.hmftools.common.purple;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.TranscriptRegion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface GeneCopyNumber extends TranscriptRegion {

    double maxCopyNumber();
    double minCopyNumber();

    int somaticRegions();
    int minRegions();

    int minRegionStart();
    int minRegionEnd();
    int depthWindowCount();

    default int minRegionBases() {
        return minRegionEnd() - minRegionStart() + 1;
    }

    SegmentSupport minRegionStartSupport();

    SegmentSupport minRegionEndSupport();

    CopyNumberMethod minRegionMethod();

    double minMinorAlleleCopyNumber();

    default int totalRegions() {
        return somaticRegions();
    }

    static Map<String,List<GeneCopyNumber>> listToMap(final List<GeneCopyNumber> geneCopyNumbers)
    {
        final Map<String,List<GeneCopyNumber>> geneCopyNumberMap = Maps.newHashMap();

        for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            List<GeneCopyNumber> geneRegions = geneCopyNumberMap.get(geneCopyNumber.geneName());

            if(geneRegions == null)
            {
                geneRegions = Lists.newArrayList();
                geneCopyNumberMap.put(geneCopyNumber.geneName(), geneRegions);
            }

            geneRegions.add(geneCopyNumber);
        }

        return geneCopyNumberMap;
    }
}
