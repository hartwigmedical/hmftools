package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.HlaCommon.MHC_CLASS_I;
import static com.hartwig.hmftools.compar.ComparTestUtil.combine;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.hla.ImmutableLilacAllele;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestLilacAlleleDataBuilder
{
    public String genes = MHC_CLASS_I;
    public String allele = "A*01:01";
    public int index = 0;
    public double missense = 0;
    public double nonsenseOrFrameshift = 0;
    public double splice = 0;
    public double inframeIndel = 0;
    public double synonymous = 0;
    public double tumorCopyNumber = 1;
    public int refTotal = 300;
    public int tumorTotal = 700;

    private static final Consumer<TestLilacAlleleDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.genes = "HLA_DPB1";
        b.allele = "A*01:02";
        b.index = 1;
        b.missense = 1;
        b.nonsenseOrFrameshift = 1;
        b.splice = 1;
        b.inframeIndel = 1;
        b.synonymous = 1;
        b.tumorCopyNumber = 2;
        b.refTotal = 250;
        b.tumorTotal = 650;
    };

    public static final TestComparableItemBuilder<TestLilacAlleleDataBuilder, LilacAlleleData> BUILDER =
            new TestComparableItemBuilder<>(TestLilacAlleleDataBuilder::new, TestLilacAlleleDataBuilder::build, ALTERNATE_INITIALIZER);

    public static LilacAlleleData buildFrom(final LilacAllele lilacAllele, final Consumer<TestLilacAlleleDataBuilder> initializer)
    {
        Consumer<TestLilacAlleleDataBuilder> copyInitializer = b -> {
            b.allele = lilacAllele.allele();
            b.missense = lilacAllele.somaticMissense();
            b.nonsenseOrFrameshift = lilacAllele.somaticNonsenseOrFrameshift();
            b.splice = lilacAllele.somaticSplice();
            b.inframeIndel = lilacAllele.somaticInframeIndel();
            b.synonymous = lilacAllele.somaticSynonymous();
            b.tumorCopyNumber = lilacAllele.tumorCopyNumber();
            b.refTotal = lilacAllele.refFragments();
            b.tumorTotal = lilacAllele.tumorFragments();
        };
        Consumer<TestLilacAlleleDataBuilder> combinedInitializer = combine(copyInitializer, initializer);
        return BUILDER.create(combinedInitializer);
    }

    private LilacAlleleData build()
    {
        LilacAllele alleleObject = ImmutableLilacAllele.builder()
                .genes(genes)
                .allele(allele)
                .refFragments(refTotal)
                .refUnique(-1)
                .refShared(-1)
                .refWild(-1)
                .tumorFragments(tumorTotal)
                .tumorUnique(-1)
                .tumorShared(-1)
                .tumorWild(-1)
                .rnaFragments(-1)
                .rnaUnique(-1)
                .rnaShared(-1)
                .rnaWild(-1)
                .tumorCopyNumber(tumorCopyNumber)
                .somaticMissense(missense)
                .somaticNonsenseOrFrameshift(nonsenseOrFrameshift)
                .somaticSplice(splice)
                .somaticSynonymous(synonymous)
                .somaticInframeIndel(inframeIndel)
                .build();
        return new LilacAlleleData(alleleObject, index);
    }
}
