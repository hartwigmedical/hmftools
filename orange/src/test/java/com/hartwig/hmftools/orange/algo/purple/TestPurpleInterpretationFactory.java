package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.datamodel.purple.*;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;

import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCharacteristics;
import org.jetbrains.annotations.NotNull;

public final class TestPurpleInterpretationFactory {

    private TestPurpleInterpretationFactory() {
    }

    @NotNull
    public static PurpleInterpretedData createMinimalTestPurpleData() {
        return builder().build();
    }

    @NotNull
    public static ImmutablePurpleInterpretedData.Builder builder() {
        return ImmutablePurpleInterpretedData.builder()
                .fit(createMinimalTestFitData())
                .characteristics(createMinimalTestCharacteristicsData());
    }

    @NotNull
    public static PurpleFit createMinimalTestFitData() {
        return ImmutablePurpleFit.builder()
                .qc(qcPass())
                .fittedPurityMethod(PurpleFittedPurityMethod.NORMAL)
                .purity(0D)
                .minPurity(0D)
                .maxPurity(0D)
                .containsTumorCells(true)
                .hasSufficientQuality(true)
                .ploidy(0D)
                .minPloidy(0D)
                .maxPloidy(0D)
                .build();
    }

    @NotNull
    private static PurpleQC qcPass() {
        return ImmutablePurpleQC.builder()
                .addStatus(PurpleQCStatus.PASS)
                .addGermlineAberrations(PurpleGermlineAberration.NONE)
                .amberMeanDepth(0)
                .unsupportedCopyNumberSegments(0)
                .deletedGenes(0)
                .contamination(0D)
                .build();
    }

    @NotNull
    private static PurpleCharacteristics createMinimalTestCharacteristicsData() {
        return ImmutablePurpleCharacteristics.builder()
                .wholeGenomeDuplication(false)
                .microsatelliteIndelsPerMb(0D)
                .microsatelliteStatus(PurpleMicrosatelliteStatus.UNKNOWN)
                .tumorMutationalBurdenPerMb(0D)
                .tumorMutationalBurdenStatus(PurpleTumorMutationalStatus.UNKNOWN)
                .tumorMutationalLoad(0)
                .tumorMutationalLoadStatus(PurpleTumorMutationalStatus.UNKNOWN)
                .svTumorMutationalBurden(0)
                .build();
    }
}
