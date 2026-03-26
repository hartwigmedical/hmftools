package com.hartwig.hmftools.orange.algo.purple;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogTestFactory;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.ChrArmCopyNumber;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.purple.MicrosatelliteStatus;
import com.hartwig.hmftools.common.segmentation.Arm;
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

        return ImmutablePurpleData.builder()
                .purityContext(minimalContext)
                .variantPlots(Lists.newArrayList());
    }

    private static FittedPurity emptyFit()
    {
        return new FittedPurity(0, 0, 0, 0, 0, 0);
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

    protected static ChrArmCopyNumber createArmCopyNumber(final String chromosome, final Arm arm)
    {
        return new ChrArmCopyNumber(
                HumanChromosome.fromString(chromosome), arm, 2, 2, 2, 2);

    }

}
