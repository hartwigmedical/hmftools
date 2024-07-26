package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.ImmutableTumorStats;
import com.hartwig.hmftools.datamodel.purple.PurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleFittedPurityMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleQC;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.TumorStats;

import org.jetbrains.annotations.NotNull;

public final class TestPurpleInterpretationFactory
{
    @NotNull
    public static PurpleRecord createMinimalTestPurpleData()
    {
        return builder().build();
    }

    @NotNull
    public static ImmutablePurpleRecord.Builder builder()
    {
        return ImmutablePurpleRecord.builder()
                .fit(createMinimalTestFitData())
                .characteristics(createMinimalTestCharacteristicsData())
                .tumorStats(createMinimalTumorStats());
    }

    @NotNull
    public static PurpleFit createMinimalTestFitData()
    {
        return ImmutablePurpleFit.builder()
                .qc(qcPass())
                .fittedPurityMethod(PurpleFittedPurityMethod.NORMAL)
                .purity(0D)
                .minPurity(0D)
                .maxPurity(0D)
                .ploidy(0D)
                .minPloidy(0D)
                .maxPloidy(0D)
                .build();
    }

    @NotNull
    private static PurpleQC qcPass()
    {
        return TestPurpleQCFactory.builder()
                .addStatus(PurpleQCStatus.PASS)
                .addGermlineAberrations(PurpleGermlineAberration.NONE)
                .build();
    }

    @NotNull
    private static PurpleCharacteristics createMinimalTestCharacteristicsData()
    {
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

    @NotNull
    private static TumorStats createMinimalTumorStats()
    {
        return ImmutableTumorStats.builder()
                .hotspotMutationCount(0)
                .hotspotStructuralVariantCount(0)
                .smallVariantCount(0)
                .structuralVariantsCount(0)
                .sumBafCounts(0)
                .build();
    }

    @NotNull
    private static TumorStats createTestTumorStats()
    {
        return ImmutableTumorStats.builder()
                .hotspotMutationCount(5)
                .hotspotStructuralVariantCount(1)
                .smallVariantCount(1500)
                .structuralVariantsCount(1100)
                .sumBafCounts(5000)
                .build();
    }
}
