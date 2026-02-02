package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.LilacAllele.MHC_CLASS_I;

import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.hla.ImmutableLilacQcData;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestLilacDataBuilder
{
    public String qcStatus = "PASS";
    public int totalFragments = 1600;
    public int fittedFragments = 1500;
    public int discardedIndels = 1;
    public int discardedAlignmentFragments = 30;
    public String hlaYAllele = "NONE";
    public List<LilacAllele> alleles = List.of(
            TestLilacAlleleBuilder.BUILDER.create(b -> b.allele = "A*01:01"),
            TestLilacAlleleBuilder.BUILDER.create(b -> b.allele = "A*01:01"),
            TestLilacAlleleBuilder.BUILDER.createWithAlternateDefaults(b -> b.allele = "B*01:01"),
            TestLilacAlleleBuilder.BUILDER.createWithAlternateDefaults(b -> b.allele = "B*01:02"),
            TestLilacAlleleBuilder.BUILDER.create(b -> b.allele = "C*02:01"),
            TestLilacAlleleBuilder.BUILDER.createWithAlternateDefaults(b -> b.allele = "C*03:04")
    );

    private static final Consumer<TestLilacDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.qcStatus = "WARN_LOW_COVERAGE";
        b.totalFragments = 1500;
        b.fittedFragments = 1400;
        b.discardedIndels = 30;
        b.discardedAlignmentFragments = 50;
        b.hlaYAllele = "Y*02:01";
        b.alleles = List.of(
                TestLilacAlleleBuilder.BUILDER.createWithAlternateDefaults(c -> c.allele = "A*01:01"),
                TestLilacAlleleBuilder.BUILDER.createWithAlternateDefaults(c -> c.allele = "A*01:01"),
                TestLilacAlleleBuilder.BUILDER.create(c -> c.allele = "B*01:01"),
                TestLilacAlleleBuilder.BUILDER.create(c -> c.allele = "B*01:02"),
                TestLilacAlleleBuilder.BUILDER.createWithAlternateDefaults(c -> c.allele = "C*02:01"),
                TestLilacAlleleBuilder.BUILDER.create(c -> c.allele = "C*04:05")
        );
    };

    public static final TestComparableItemBuilder<TestLilacDataBuilder, LilacData> BUILDER =
            new TestComparableItemBuilder<>(TestLilacDataBuilder::new, TestLilacDataBuilder::build, ALTERNATE_INITIALIZER);

    private LilacData build()
    {
        LilacQcData qcData = ImmutableLilacQcData.builder()
                .genes(MHC_CLASS_I)
                .status(qcStatus)
                .totalFragments(totalFragments)
                .fittedFragments(fittedFragments)
                .discardedIndels(discardedIndels)
                .discardedAlignmentFragments(discardedAlignmentFragments)
                .hlaYAllele(hlaYAllele)
                .build();
        return new LilacData(qcData, alleles);
    }
}
