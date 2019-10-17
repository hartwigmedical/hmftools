package com.hartwig.hmftools.patientreporter.structural;

import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableReportableGeneFusion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class SvAnalysisDatamodelTestFactory {

    private SvAnalysisDatamodelTestFactory() {
    }

    @NotNull
    static ImmutableReportableGeneFusion.Builder createTestFusionBuilder() {
        return ImmutableReportableGeneFusion.builder()
                .geneContextStart(Strings.EMPTY)
                .geneContextEnd(Strings.EMPTY)
                .geneTranscriptStart(Strings.EMPTY)
                .geneTranscriptEnd(Strings.EMPTY)
                .ploidy(1D);
    }

    @NotNull
    static ImmutableReportableDisruption.Builder createTestDisruptionBuilder() {
        return ImmutableReportableDisruption.builder()
                .svId(0)
                .chromosome(Strings.EMPTY)
                .orientation(1)
                .type(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .chrBand(Strings.EMPTY)
                .strand(1)
                .exonUp(0)
                .exonDown(0)
                .undisruptedCopyNumber(1);
    }
}
