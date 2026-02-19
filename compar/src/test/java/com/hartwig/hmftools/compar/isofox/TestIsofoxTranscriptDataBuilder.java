package com.hartwig.hmftools.compar.isofox;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.rna.ImmutableTranscriptExpression;
import com.hartwig.hmftools.common.rna.TranscriptExpression;
import com.hartwig.hmftools.compar.TestComparableItemBuilder;

public class TestIsofoxTranscriptDataBuilder
{
    public String transcriptName = "ENST00000304494";
    public String geneName = "CDKN2A";
    public double tpm = 2.0;

    private static final Consumer<TestIsofoxTranscriptDataBuilder> ALTERNATE_INITIALIZER = b ->
    {
        b.transcriptName = "ENST00000646891";
        b.geneName = "BRAF";
        b.tpm = 0.001;
    };

    public static final TestComparableItemBuilder<TestIsofoxTranscriptDataBuilder, IsofoxTranscriptData> BUILDER =
            new TestComparableItemBuilder<>(TestIsofoxTranscriptDataBuilder::new, TestIsofoxTranscriptDataBuilder::build, ALTERNATE_INITIALIZER);

    private IsofoxTranscriptData build()
    {
        final TranscriptExpression transcriptExpression = ImmutableTranscriptExpression.builder()
                .transcriptName(transcriptName)
                .geneName(geneName)
                .tpm(tpm)
                .build();
        return new IsofoxTranscriptData(transcriptExpression);
    }
}
