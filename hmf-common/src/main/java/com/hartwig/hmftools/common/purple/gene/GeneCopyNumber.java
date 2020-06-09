package com.hartwig.hmftools.common.purple.gene;

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

    SegmentSupport minRegionStartSupport();

    SegmentSupport minRegionEndSupport();

    CopyNumberMethod minRegionMethod();

    double minMinorAlleleCopyNumber();

    default int totalRegions() {
        return somaticRegions() + germlineHet2HomRegions() + germlineHomRegions();
    }
}
