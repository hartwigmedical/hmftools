package com.hartwig.hmftools.compar.isofox;

import java.util.Collections;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.rna.ImmutableTranscriptExpression;
import com.hartwig.hmftools.common.rna.TranscriptExpression;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestRnaTranscriptDataBuilder
{
    public String transcriptName = "ENST00000304494";
    public String geneName = "CDKN2A";
    public double tpm = 2.0;

    private static final Consumer<TestRnaTranscriptDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.transcriptName = "ENST00000646891";
        b.geneName = "BRAF";
        b.tpm = 0.001;
    };

    public static final TestComparableItemBuilder<TestRnaTranscriptDataBuilder, RnaTranscriptData> BUILDER =
            new TestComparableItemBuilder<>(TestRnaTranscriptDataBuilder::new, TestRnaTranscriptDataBuilder::build, ALTERNATE_INITIALIZER);

    private RnaTranscriptData build()
    {
        final TranscriptExpression transcriptExpression = ImmutableTranscriptExpression.builder()
                .transcriptName(transcriptName)
                .geneName(geneName)
                .tpm(tpm)
                .build();
        return new RnaTranscriptData(
                transcriptExpression, new RnaTranscriptDataComparer(null, Collections.emptyMap()).fieldsList());
    }
}
