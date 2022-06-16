package com.hartwig.hmftools.common.linx;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.sv.linx.FusionPhasedType;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class LinxTestFactory {

    private LinxTestFactory() {
    }

    @NotNull
    public static LinxFusion createMinimalTestFusion() {
        return fusionBuilder().build();
    }

    @NotNull
    public static ImmutableLinxFusion.Builder fusionBuilder() {
        return ImmutableLinxFusion.builder()
                .fivePrimeBreakendId(0)
                .threePrimeBreakendId(0)
                .name(Strings.EMPTY)
                .reported(false)
                .reportedType(KnownFusionType.NONE.toString())
                .phased(FusionPhasedType.OUT_OF_FRAME)
                .likelihood(FusionLikelihoodType.NA)
                .chainLength(0)
                .chainLinks(0)
                .chainTerminated(false)
                .domainsKept(Strings.EMPTY)
                .domainsLost(Strings.EMPTY)
                .skippedExonsUp(0)
                .skippedExonsDown(0)
                .fusedExonUp(0)
                .fusedExonDown(0)
                .geneStart(Strings.EMPTY)
                .geneContextStart(Strings.EMPTY)
                .geneTranscriptStart(Strings.EMPTY)
                .geneEnd(Strings.EMPTY)
                .geneContextEnd(Strings.EMPTY)
                .geneTranscriptEnd(Strings.EMPTY)
                .junctionCopyNumber(0D);
    }

    @NotNull
    public static ImmutableLinxBreakend.Builder breakendBuilder() {
        return ImmutableLinxBreakend.builder()
                .id(0)
                .svId(0)
                .isStart(true)
                .gene(Strings.EMPTY)
                .transcriptId(Strings.EMPTY)
                .canonical(true)
                .geneOrientation(Strings.EMPTY)
                .disruptive(false)
                .reportedDisruption(false)
                .undisruptedCopyNumber(0D)
                .regionType(Strings.EMPTY)
                .codingContext(Strings.EMPTY)
                .biotype(Strings.EMPTY)
                .exonicBasePhase(0)
                .nextSpliceExonRank(0)
                .nextSpliceExonPhase(0)
                .nextSpliceDistance(0)
                .totalExonCount(0)
                .type(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .orientation(0)
                .strand(0)
                .chrBand(Strings.EMPTY)
                .exonUp(0)
                .exonDown(0)
                .junctionCopyNumber(0D);
    }
}
