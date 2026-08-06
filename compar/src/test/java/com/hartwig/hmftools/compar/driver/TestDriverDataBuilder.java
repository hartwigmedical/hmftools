package com.hartwig.hmftools.compar.driver;

import java.util.Collections;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.MicrosatelliteStatus;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestDriverDataBuilder
{
    public String gene = "CDKN2A";
    public DriverType driver = DriverType.DEL;
    public String transcript = "ENST00000579755";
    public boolean isCanonical = false;
    public LikelihoodMethod likelihoodMethod = LikelihoodMethod.DEL;
    public double likelihood = 0.6;
    public double minCopyNumber = 0.1;
    public double maxCopyNumber = 0.1;
    public String chromosomeBand = "9p21";
    public String comparisonChromosome = "chr9";
    public boolean checkTranscript = true;
    public boolean isPass = true;

    private static final Consumer<TestDriverDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.gene = "BRAF";
        b.driver = DriverType.MUTATION;
        b.transcript = "ENST00000646891";
        b.isCanonical = true;
        b.likelihoodMethod = LikelihoodMethod.HOTSPOT;
        b.likelihood = 0.9;
        b.minCopyNumber = 0.6;
        b.maxCopyNumber = 0.9;
        b.chromosomeBand = "7q34";
        b.comparisonChromosome = "chr7";
        b.checkTranscript = false;
        b.isPass = true;
    };

    public static final TestComparableItemBuilder<TestDriverDataBuilder, DriverData> BUILDER =
            new TestComparableItemBuilder<>(TestDriverDataBuilder::new, TestDriverDataBuilder::build, ALTERNATE_INITIALIZER);

    private DriverData build()
    {
        final DriverCatalog driverCatalog = ImmutableDriverCatalog.builder()
                .chromosomeBand(chromosomeBand)
                .gene(gene)
                .transcript(transcript)
                .isCanonical(isCanonical)
                .driver(driver)
                .likelihoodMethod(likelihoodMethod)
                .driverLikelihood(likelihood)
                .reportedStatus(ReportedStatus.REPORTED)
                .minCopyNumber(minCopyNumber)
                .maxCopyNumber(maxCopyNumber)
                .category(DriverCategory.ONCO)
                .chromosome("")
                .missense(-1)
                .nonsense(-1)
                .splice(-1)
                .inframe(-1)
                .frameshift(-1)
                .biallelic(true)
                .build();

        return new DriverData(
                driverCatalog, buildPurityData(), comparisonChromosome, checkTranscript, isPass,
                new DriverComparer(null, Collections.emptyMap()).fieldsList());
    }

    protected static PurplePurity buildPurityData()
    {
        return new PurplePurity(
                0.5, 1, 1, 1, 2, Gender.FEMALE, FittedPurityMethod.NORMAL, 1,
                1, 1, 1, 1, 1, 1,
                0, false, 0, MicrosatelliteStatus.MSS, 0, TumorMutationalStatus.LOW,
                0, TumorMutationalStatus.LOW, 0, RunMode.TUMOR_GERMLINE, false);
    }

    protected static PurplePurity buildAlternatePurityData()
    {
        return new PurplePurity(
                0.8, 1, 1, 1, 3, Gender.FEMALE, FittedPurityMethod.NORMAL, 1,
                1, 1, 1, 1, 1, 1,
                0, false, 0, MicrosatelliteStatus.MSS, 0, TumorMutationalStatus.LOW,
                0, TumorMutationalStatus.LOW, 0, RunMode.TUMOR_GERMLINE, false);
    }
}
