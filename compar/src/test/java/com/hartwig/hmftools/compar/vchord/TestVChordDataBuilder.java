package com.hartwig.hmftools.compar.vchord;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.vchord.ImmutableVChordPrediction;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestVChordDataBuilder
{
    public double breast = 0.1;
    public double ovarian = 0.05;
    public double pancreatic = 0.15;
    public double prostate = 0.2;
    public double other = 0.25;

    private static final Consumer<TestVChordDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.breast = 0.4;
        b.ovarian = 0.6;
        b.pancreatic = 0.9;
        b.prostate = 0.85;
        b.other = 1.0;
    };

    public static final TestComparableItemBuilder<TestVChordDataBuilder, VChordData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestVChordDataBuilder::new,
                    TestVChordDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private VChordData build()
    {
        return new VChordData(ImmutableVChordPrediction.builder()
                .breastCancerHrdScore(breast)
                .ovarianCancerHrdScore(ovarian)
                .pancreaticCancerScore(pancreatic)
                .prostateCancerScore(prostate)
                .otherCancerScore(other)
                .build());
    }
}
