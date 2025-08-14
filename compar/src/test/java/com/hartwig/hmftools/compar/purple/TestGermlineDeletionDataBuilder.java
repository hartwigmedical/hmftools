package com.hartwig.hmftools.compar.purple;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDetectionMethod;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestGermlineDeletionDataBuilder
{
    public String gene = "BRAF";
    public boolean reported = true;
    public GermlineStatus germlineStatus = GermlineStatus.HET_DELETION;
    public GermlineStatus tumorStatus = GermlineStatus.HOM_DELETION;
    public double germlineCopyNumber = 1;
    public double tumorCopyNumber = 0;
    public String comparisonChromosome = "chr7";
    public String chromosomeBand = "7q34";

    private static final Consumer<TestGermlineDeletionDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.gene = "BRCA2";
        b.reported = false;
        b.germlineStatus = GermlineStatus.DIPLOID;
        b.tumorStatus = GermlineStatus.AMPLIFICATION;
        b.germlineCopyNumber = 2;
        b.tumorCopyNumber = 3;
        b.comparisonChromosome = "chr13";
        b.chromosomeBand = "13q13.1";
    };

    public static final TestComparableItemBuilder<TestGermlineDeletionDataBuilder, GermlineDeletionData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestGermlineDeletionDataBuilder::new,
                    TestGermlineDeletionDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private GermlineDeletionData build()
    {
        return new GermlineDeletionData(
                new GermlineDeletion(
                        gene,
                        "",
                        chromosomeBand,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        GermlineDetectionMethod.SEGMENT,
                        germlineStatus,
                        tumorStatus,
                        germlineCopyNumber,
                        tumorCopyNumber,
                        "",
                        -1,
                        reported
                ),
                comparisonChromosome
        );
    }
}
