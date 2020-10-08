package com.hartwig.hmftools.protect.structural;

import com.hartwig.hmftools.common.variant.structural.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.variant.structural.linx.FusionPhasedType;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxFusion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

final class SvAnalysisDatamodelTestFactory {

    private SvAnalysisDatamodelTestFactory() {
    }

    @NotNull
    static ImmutableLinxFusion.Builder createTestFusionBuilder() {
        return ImmutableLinxFusion.builder()
                .fivePrimeBreakendId(1)
                .threePrimeBreakendId(2)
                .name(Strings.EMPTY)
                .reported(true)
                .reportedType(Strings.EMPTY)
                .phased(FusionPhasedType.SKIPPED_EXONS)
                .likelihood(FusionLikelihoodType.HIGH)
                .chainLength(1)
                .chainLinks(1)
                .chainTerminated(true)
                .domainsKept(Strings.EMPTY)
                .domainsLost(Strings.EMPTY)
                .skippedExonsUp(2)
                .skippedExonsDown(4)
                .fusedExonUp(6)
                .fusedExonDown(7)
                .geneStart(Strings.EMPTY)
                .geneContextStart(Strings.EMPTY)
                .geneTranscriptStart(Strings.EMPTY)
                .geneEnd(Strings.EMPTY)
                .geneContextEnd(Strings.EMPTY)
                .geneTranscriptEnd(Strings.EMPTY)
                .junctionCopyNumber(1D);
    }

    @NotNull
    static ImmutableLinxBreakend.Builder createTestDisruptionBuilder() {
        return ImmutableLinxBreakend.builder()
                .id(0)
                .svId(0)
                .isStart(true)
                .gene(Strings.EMPTY)
                .transcriptId(Strings.EMPTY)
                .canonical(true)
                .geneOrientation(Strings.EMPTY)
                .disruptive(true)
                .reportedDisruption(true)
                .undisruptedCopyNumber(0.1)
                .regionType(Strings.EMPTY)
                .codingContext(Strings.EMPTY)
                .biotype(Strings.EMPTY)
                .exonicBasePhase(1)
                .nextSpliceExonRank(1)
                .nextSpliceExonPhase(1)
                .nextSpliceDistance(1)
                .totalExonCount(1)
                .type(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .orientation(1)
                .strand(1)
                .chrBand(Strings.EMPTY)
                .exonUp(0)
                .exonDown(0)
                .junctionCopyNumber(0.1);
    }
}
