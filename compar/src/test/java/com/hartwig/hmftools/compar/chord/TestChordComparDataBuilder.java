package com.hartwig.hmftools.compar.chord;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.chord.ChordStatus;
import com.hartwig.hmftools.common.chord.ImmutableChordData;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestChordComparDataBuilder
{
    public double brca1 = 0.1;
    public double brca2 = 0.05;
    public double score = 0.15;
    public String type = "none";
    public ChordStatus status = ChordStatus.HR_PROFICIENT;

    private static final Consumer<TestChordComparDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.brca1 = 0.3;
        b.brca2 = 0.6 ;
        b.score = 0.9;
        b.type = "BRCA2_type";
        b.status = ChordStatus.HR_DEFICIENT;
    };

    public static final TestComparableItemBuilder<TestChordComparDataBuilder, ChordComparData> BUILDER =
            new TestComparableItemBuilder<>(
                    TestChordComparDataBuilder::new,
                    TestChordComparDataBuilder::build,
                    ALTERNATE_INITIALIZER
            );

    private ChordComparData build()
    {
        return new ChordComparData(ImmutableChordData.builder()
                .BRCA1Value(brca1)
                .BRCA2Value(brca2)
                .hrdValue(score)
                .hrdType(type)
                .hrStatus(status)
                .remarksHrStatus("")
                .remarksHrdType("")
                .build());
    }
}
