package com.hartwig.hmftools.protect;

import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class PurpleTestFactory {

    private PurpleTestFactory() {
    }

    @NotNull
    public static ReportableGainLoss testReportableGainLoss(@NotNull String gene, @NotNull CopyNumberInterpretation interpretation) {
        return ImmutableReportableGainLoss.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .interpretation(interpretation)
                .copies(1)
                .build();
    }
}
