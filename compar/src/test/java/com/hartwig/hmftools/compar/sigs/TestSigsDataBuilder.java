package com.hartwig.hmftools.compar.sigs;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.sigs.ImmutableSignatureAllocation;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestSigsDataBuilder
{
    public String signature = "Sig1";
    public double percent = 0.1;

    private static final Consumer<TestSigsDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.signature = "Sig2";
        b.percent = 0.9;
    };

    public static final TestComparableItemBuilder<TestSigsDataBuilder, SigsData> BUILDER =
            new TestComparableItemBuilder<>(TestSigsDataBuilder::new, TestSigsDataBuilder::build, ALTERNATE_INITIALIZER);

    private SigsData build()
    {
        final SignatureAllocation signatureAllocation = ImmutableSignatureAllocation.builder()
                .signature(signature)
                .allocation(-1)
                .percent(percent)
                .build();
        return new SigsData(signatureAllocation);
    }
}
