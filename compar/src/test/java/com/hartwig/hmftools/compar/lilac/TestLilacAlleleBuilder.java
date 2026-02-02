package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.HlaCommon.MHC_CLASS_I;
import static com.hartwig.hmftools.compar.ComparTestUtil.combine;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.hla.ImmutableLilacAllele;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestLilacAlleleBuilder
{
    public String allele = "A*01:01";
    public double missense = 0;
    public double nonsenseOrFrameshift = 0;
    public double splice = 0;
    public double inframeIndel = 0;
    public double synonymous = 0;
    public double tumorCopyNumber = 1;
    public int refTotal = 300;
    public int tumorTotal = 700;

    private static final Consumer<TestLilacAlleleBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.allele = "A*01:02";
        b.missense = 1;
        b.nonsenseOrFrameshift = 1;
        b.splice = 1;
        b.inframeIndel = 1;
        b.synonymous = 1;
        b.tumorCopyNumber = 2;
        b.refTotal = 250;
        b.tumorTotal = 650;
    };

    public static final TestComparableItemBuilder<TestLilacAlleleBuilder, LilacAllele> BUILDER =
            new TestComparableItemBuilder<>(TestLilacAlleleBuilder::new, TestLilacAlleleBuilder::build, ALTERNATE_INITIALIZER);

    public static LilacAllele buildFrom(final LilacAllele lilacAllele, final Consumer<TestLilacAlleleBuilder> initializer)
    {
        Consumer<TestLilacAlleleBuilder> copyInitializer = b -> {
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
        Consumer<TestLilacAlleleBuilder> combinedInitializer = combine(copyInitializer, initializer);
        return BUILDER.create(combinedInitializer);
    }

    private LilacAllele build()
    {
        return ImmutableLilacAllele.builder()
                .genes(MHC_CLASS_I)
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
    }
}
