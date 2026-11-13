package com.hartwig.hmftools.compar.lilac;

import static com.hartwig.hmftools.common.hla.HlaCommon.MHC_CLASS_I;

import java.util.Collections;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.hla.ImmutableLilacQcData;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestLilacQcDataBuilder
{
    public String genes = MHC_CLASS_I;
    public String qcStatus = "PASS";
    public int totalFragments = 1600;
    public int fittedFragments = 1500;
    public int discardedIndels = 1;
    public int discardedAlignmentFragments = 30;
    public String hlaYAllele = "NONE";

    public String allelesStr = "A*01:01;A*01:01;B*01:01;B*01:02;C*02:01;C*03:04";

    private static final Consumer<TestLilacQcDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.genes = "HLA_DPB1";
        b.qcStatus = "WARN_LOW_COVERAGE";
        b.totalFragments = 1500;
        b.fittedFragments = 1400;
        b.discardedIndels = 30;
        b.discardedAlignmentFragments = 50;
        b.hlaYAllele = "Y*02:01";

        b.allelesStr = "A*01:01;A*01:01;B*01:01;B*01:02;C*02:01;C*04:05";
    };

    public static final TestComparableItemBuilder<TestLilacQcDataBuilder, LilacQcComparData> BUILDER =
            new TestComparableItemBuilder<>(TestLilacQcDataBuilder::new, TestLilacQcDataBuilder::build, ALTERNATE_INITIALIZER);

    private LilacQcComparData build()
    {
        com.hartwig.hmftools.common.hla.LilacQcData qcData = ImmutableLilacQcData.builder()
                .genes(genes)
                .status(qcStatus)
                .totalFragments(totalFragments)
                .fittedFragments(fittedFragments)
                .discardedIndels(discardedIndels)
                .discardedAlignmentFragments(discardedAlignmentFragments)
                .hlaYAllele(hlaYAllele)
                .build();

        return new LilacQcComparData(qcData, allelesStr, new LilacQcComparer(null, Collections.emptyMap()).fieldsList());
    }
}
