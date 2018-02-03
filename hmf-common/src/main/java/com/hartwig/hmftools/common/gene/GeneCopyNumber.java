package com.hartwig.hmftools.common.gene;

import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface GeneCopyNumber extends TranscriptRegion, CopyNumber {

    double maxCopyNumber();

    double minCopyNumber();

    double meanCopyNumber();

    int somaticRegions();

    int germlineHet2HomRegions();

    int germlineHomRegions();

    int minRegions();

    long minRegionStart();

    long minRegionEnd();

    SegmentSupport minRegionStartSupport();

    SegmentSupport minRegionEndSupport();

    CopyNumberMethod minRegionMethod();

    int nonsenseBiallelicCount();

    int nonsenseNonBiallelicCount();

    double nonsenseNonBiallelicPloidy();

    int spliceBiallelicCount();

    int spliceNonBiallelicCount();

    double spliceNonBiallelicPloidy();

    int missenseBiallelicCount();

    int missenseNonBiallelicCount();

    double missenseNonBiallelicPloidy();

    default int totalRegions() {
        return somaticRegions() + germlineHet2HomRegions() + germlineHomRegions();
    }

    @Override
    default int value() {
        return (int) Math.max(0, Math.round(minCopyNumber()));
    }
}
