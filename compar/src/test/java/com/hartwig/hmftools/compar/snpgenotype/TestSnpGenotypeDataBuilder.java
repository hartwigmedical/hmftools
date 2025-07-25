package com.hartwig.hmftools.compar.snpgenotype;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestSnpGenotypeDataBuilder
{
    public String chromosome = "7";
    public int position = 10000;
    public String ref = "A";
    public String alt = "T";
    public String genotype = "HET";
    public String vcfSampleId = "FAKE_ID";
    public String comparisonChromosome = "7";
    public int comparisonPosition = 10000;

    private static final Consumer<TestSnpGenotypeDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.chromosome = "10";
        b.position = 20000;
        b.ref = "C";
        b.alt = "G";
        b.genotype = "HOM_VAR";
        b.vcfSampleId = "NA";
        b.comparisonChromosome = "10";
        b.comparisonPosition = 20000;
    };

    public static final TestComparableItemBuilder<TestSnpGenotypeDataBuilder, SnpGenotypeData> BUILDER =
            new TestComparableItemBuilder<>(TestSnpGenotypeDataBuilder::new, TestSnpGenotypeDataBuilder::build, ALTERNATE_INITIALIZER);

    private SnpGenotypeData build()
    {
        BasePosition comparisonBasePosition = new BasePosition(comparisonChromosome, comparisonPosition);
        return new SnpGenotypeData(chromosome, position, ref, alt, genotype, vcfSampleId, comparisonBasePosition);
    }
}
