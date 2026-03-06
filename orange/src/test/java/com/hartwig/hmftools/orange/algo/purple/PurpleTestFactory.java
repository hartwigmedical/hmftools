package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogTestFactory;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.purple.MicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

public final class PurpleTestFactory
{
    public static PurpleData createMinimalTestPurpleData()
    {
        return createMinimalTestPurpleDataBuilder().build();
    }

    public static ImmutablePurpleDriver.Builder purpleDriverBuilder()
    {
        DriverCatalog catalog = DriverCatalogTestFactory.builder().build();
        PurpleDriver driver = PurpleConversion.convert(catalog);
        return ImmutablePurpleDriver.builder().from(driver);
    }

    public static ImmutablePurpleData.Builder createMinimalTestPurpleDataBuilder()
    {
        PurityContext minimalContext = ImmutablePurityContext.builder()
                .gender(Gender.FEMALE)
                .runMode(RunMode.TUMOR_GERMLINE)
                .targeted(false)
                .bestFit(emptyFit())
                .method(FittedPurityMethod.NORMAL)
                .score(emptyScore())
                .qc(qcPass())
                .polyClonalProportion(0D)
                .wholeGenomeDuplication(false)
                .microsatelliteIndelsPerMb(0D)
                .tumorMutationalBurdenPerMb(0D)
                .tumorMutationalLoad(0)
                .svTumorMutationalBurden(0)
                .microsatelliteStatus(MicrosatelliteStatus.UNKNOWN)
                .tumorMutationalLoadStatus(TumorMutationalStatus.UNKNOWN)
                .tumorMutationalBurdenStatus(TumorMutationalStatus.UNKNOWN)
                .build();

        return ImmutablePurpleData.builder().purityContext(minimalContext);
    }

    private static FittedPurity emptyFit()
    {
        return ImmutableFittedPurity.builder()
                .purity(0D)
                .normFactor(0D)
                .ploidy(0D)
                .score(0D)
                .diploidProportion(0D)
                .somaticPenalty(0D)
                .build();
    }

    private static FittedPurityScore emptyScore()
    {
        return ImmutableFittedPurityScore.builder()
                .minPurity(0D)
                .maxPurity(0D)
                .minPloidy(0D)
                .maxPloidy(0D)
                .minDiploidProportion(0D)
                .maxDiploidProportion(0D)
                .build();
    }

    private static PurpleQC qcPass()
    {
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
                .lohPercent(0)
                .tincLevel(0)
                .chimerismPercentage(0)
                .build();
    }
}
