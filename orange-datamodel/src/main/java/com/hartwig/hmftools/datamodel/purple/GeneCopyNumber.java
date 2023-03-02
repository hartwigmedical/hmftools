package com.hartwig.hmftools.datamodel.purple;

import com.hartwig.hmftools.datamodel.genome.region.TranscriptRegion;
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
}
