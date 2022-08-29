package com.hartwig.hmftools.common.purple;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class GeneCopyNumberTestFactory {

    private GeneCopyNumberTestFactory() {
    }

    @NotNull
    public static ImmutableGeneCopyNumber.Builder builder() {
        return ImmutableGeneCopyNumber.builder()
                .start(0)
                .end(0)
                .geneName(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(1)
                .somaticRegions(1)
                .minCopyNumber(0D)
                .maxCopyNumber(0D)
                .transName(Strings.EMPTY)
                .isCanonical(true)
                .minMinorAlleleCopyNumber(0)
                .bafWindows(0);
    }
}
