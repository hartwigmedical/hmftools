package com.hartwig.hmftools.compar.peach;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.peach.ImmutablePeachGenotype;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestPeachDataBuilder
{
    public String gene = "DPYD";
    public String allele = "*2A";
    public int alleleCount = 2;
    public String function = "No Function";
    public String linkedDrugs = "5-Fluorouracil;Capecitabine;Tegafur";
    public String prescriptionUrls = "https://www.pharmgkb.org/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/guidelineAnnotation/PA166104963;https://www.pharmgkb.org/guidelineAnnotation/PA166104944";

    private static final Consumer<TestPeachDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.gene = "UGT1A1";
        b.allele = "*28";
        b.alleleCount = 1;
        b.function = "Reduced Function";
        b.linkedDrugs = "Irinotecan";
        b.prescriptionUrls = "https://www.pharmgkb.org/guidelineAnnotation/PA166104951";
    };

    public static final TestComparableItemBuilder<TestPeachDataBuilder, PeachData> BUILDER =
            new TestComparableItemBuilder<>(TestPeachDataBuilder::new, TestPeachDataBuilder::build, ALTERNATE_INITIALIZER);

    private PeachData build()
    {
        return new PeachData(ImmutablePeachGenotype.builder()
                .gene(gene)
                .allele(allele)
                .alleleCount(alleleCount)
                .function(function)
                .linkedDrugs(linkedDrugs)
                .urlPrescriptionInfo(prescriptionUrls)
                .build());
    }
}
