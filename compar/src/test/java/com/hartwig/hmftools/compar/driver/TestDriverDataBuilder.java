package com.hartwig.hmftools.compar.driver;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.purple.ReportedStatus;
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
    public String chromosome = "chr9";
    public String chromosomeBand = "9p21";
    public String comparisonChromosome = "chr9";
    public boolean checkTranscript = true;

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
        b.chromosome = "chr7";
        b.chromosomeBand = "7q34";
        b.comparisonChromosome = "chr7";
        b.checkTranscript = false;
    };

    public static final TestComparableItemBuilder<TestDriverDataBuilder, DriverData> BUILDER =
            new TestComparableItemBuilder<>(TestDriverDataBuilder::new, TestDriverDataBuilder::build, ALTERNATE_INITIALIZER);

    private DriverData build()
    {
        final DriverCatalog driverCatalog = ImmutableDriverCatalog.builder()
                .chromosome(chromosome)
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
                .missense(-1)
                .nonsense(-1)
                .splice(-1)
                .inframe(-1)
                .frameshift(-1)
                .biallelic(true)
                .build();
        return new DriverData(driverCatalog, comparisonChromosome, checkTranscript);
    }
}
