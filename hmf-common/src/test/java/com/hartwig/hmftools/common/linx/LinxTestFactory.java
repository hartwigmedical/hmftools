package com.hartwig.hmftools.common.linx;

import com.hartwig.hmftools.common.sv.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.sv.linx.FusionPhasedType;
import com.hartwig.hmftools.common.sv.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class LinxTestFactory {

    private LinxTestFactory() {
    }

    @NotNull
    public static LinxFusion createMinimalTestFusion() {
        return ImmutableLinxFusion.builder()
                .fivePrimeBreakendId(0)
                .threePrimeBreakendId(0)
                .name(Strings.EMPTY)
                .reported(false)
                .reportedType(Strings.EMPTY)
                .phased(FusionPhasedType.INFRAME)
                .likelihood(FusionLikelihoodType.HIGH)
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
                .junctionCopyNumber(0D)
                .build();
    }
}
