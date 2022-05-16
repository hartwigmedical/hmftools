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
    public static PurpleData createMinimalTestPurpleData() {
        return ImmutablePurpleData.builder()
                .qc(qcPass())
                .fittedPurityMethod(FittedPurityMethod.NORMAL)
                .purity(0D)
                .minPurity(0D)
                .maxPurity(0D)
                .hasReliablePurity(true)
                .hasReliableQuality(true)
                .ploidy(0D)
                .minPloidy(0D)
                .maxPloidy(0D)
                .wholeGenomeDuplication(false)
                .microsatelliteIndelsPerMb(0D)
                .microsatelliteStatus(MicrosatelliteStatus.UNKNOWN)
                .tumorMutationalBurdenPerMb(0D)
                .tumorMutationalLoad(0)
                .tumorMutationalLoadStatus(TumorMutationalStatus.UNKNOWN)
                .svTumorMutationalBurden(0)
                .build();
    }

    @NotNull
    public static ReportableGainLoss createTestReportableGainLoss() {
        return createReportableGainLoss(Strings.EMPTY, CopyNumberInterpretation.FULL_GAIN);
    }

    @NotNull
    public static ReportableGainLoss createReportableGainLoss(@NotNull String gene, @NotNull CopyNumberInterpretation interpretation) {
        return ImmutableReportableGainLoss.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript("123")
                .isCanonical(true)
                .interpretation(interpretation)
                .minCopies(1)
                .maxCopies(1)
                .build();
    }

    @NotNull
    private static PurpleQC qcPass() {
        return ImmutablePurpleQC.builder()
                .method(FittedPurityMethod.NORMAL)
                .amberMeanDepth(0)
                .copyNumberSegments(1)
                .unsupportedCopyNumberSegments(0)
                .deletedGenes(0)
                .purity(0.5)
                .contamination(0D)
                .cobaltGender(Gender.FEMALE)
                .amberGender(Gender.FEMALE)
                .build();
    }
}
