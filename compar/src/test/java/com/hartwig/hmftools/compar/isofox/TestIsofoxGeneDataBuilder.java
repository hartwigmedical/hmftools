package com.hartwig.hmftools.compar.isofox;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.ImmutableGeneExpression;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestIsofoxGeneDataBuilder
{
    public String geneName = "CDKN2A";
    public int splicedFragments = 500;
    public int unsplicedFragments = 100;
    public double tpm = 2.0;

    private static final Consumer<TestIsofoxGeneDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.geneName = "BRAF";
        b.splicedFragments = 50;
        b.unsplicedFragments = 10;
        b.tpm = 0.001;
    };

    public static final TestComparableItemBuilder<TestIsofoxGeneDataBuilder, IsofoxGeneData> BUILDER =
            new TestComparableItemBuilder<>(TestIsofoxGeneDataBuilder::new, TestIsofoxGeneDataBuilder::build, ALTERNATE_INITIALIZER);

    private IsofoxGeneData build()
    {
        final GeneExpression geneExpression = ImmutableGeneExpression.builder()
                .geneName(geneName)
                .tpm(tpm)
                .splicedFragments(splicedFragments)
                .unsplicedFragments(unsplicedFragments)
                .medianTpmCancer(0.1)
                .percentileCancer(0.5)
                .medianTpmCohort(0.3)
                .percentileCohort(0.4)
                .build();
        return new IsofoxGeneData(geneExpression);
    }
}
