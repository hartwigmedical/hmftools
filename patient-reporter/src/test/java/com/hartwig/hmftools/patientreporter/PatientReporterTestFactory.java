package com.hartwig.hmftools.patientreporter;

import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PatientReporterTestFactory {

    private PatientReporterTestFactory() {
    }

    @NotNull
    public static ImmutableGeneCopyNumber.Builder createTestCopyNumberBuilder() {
        return ImmutableGeneCopyNumber.builder()
                .start(1)
                .end(2)
                .gene(Strings.EMPTY)
                .chromosome("1")
                .chromosomeBand(Strings.EMPTY)
                .minRegionStart(0)
                .minRegionStartSupport(SegmentSupport.NONE)
                .minRegionEnd(0)
                .minRegionEndSupport(SegmentSupport.NONE)
                .minRegionMethod(CopyNumberMethod.UNKNOWN)
                .minRegions(1)
                .germlineHet2HomRegions(0)
                .germlineHomRegions(0)
                .somaticRegions(1)
                .minCopyNumber(0)
                .maxCopyNumber(0)
                .transcriptID("trans")
                .transcriptVersion(0)
                .minMinorAllelePloidy(0);
    }
}
