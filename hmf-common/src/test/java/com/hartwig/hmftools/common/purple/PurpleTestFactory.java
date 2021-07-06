package com.hartwig.hmftools.common.purple;

import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.copynumber.ImmutableReportableGainLoss;
import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.common.purple.purity.FittedPurityMethod;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;

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

    @NotNull
    public static ImmutablePurpleData.Builder testPurpleDataBuilder() {
        return ImmutablePurpleData.builder()
                .addPurpleQC(PurpleQCStatus.PASS)
                .fittedPurityMethod(FittedPurityMethod.NORMAL)
                .wholeGenomeDuplication(false)
                .purity(0D)
                .minPurity(0D)
                .maxPurity(0D)
                .hasReliablePurity(true)
                .hasReliableQuality(true)
                .ploidy(0D)
                .minPloidy(0D)
                .maxPloidy(0D)
                .microsatelliteIndelsPerMb(0D)
                .microsatelliteStatus(MicrosatelliteStatus.UNKNOWN)
                .tumorMutationalBurdenPerMb(0D)
                .tumorMutationalLoad(0)
                .tumorMutationalLoadStatus(TumorMutationalStatus.UNKNOWN);
    }
}
