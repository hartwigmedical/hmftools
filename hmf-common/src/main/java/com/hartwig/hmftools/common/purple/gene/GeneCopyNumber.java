package com.hartwig.hmftools.common.purple.gene;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.TranscriptRegion;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface GeneCopyNumber extends TranscriptRegion {

    double maxCopyNumber();

    double minCopyNumber();

    int somaticRegions();

    int germlineHet2HomRegions();

    int germlineHomRegions();

    int minRegions();

    long minRegionStart();

    long minRegionEnd();

    default long minRegionBases() {
        return minRegionEnd() - minRegionStart() + 1;
    }

    SegmentSupport minRegionStartSupport();

    SegmentSupport minRegionEndSupport();

    CopyNumberMethod minRegionMethod();

    double minMinorAlleleCopyNumber();

    default int totalRegions() {
        return somaticRegions() + germlineHet2HomRegions() + germlineHomRegions();
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
