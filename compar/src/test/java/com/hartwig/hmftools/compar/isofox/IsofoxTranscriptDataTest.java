package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.GeneExpressionFile.FLD_ADJ_TPM;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;
import org.junit.Test;

public class IsofoxTranscriptDataTest extends ComparableItemTest<IsofoxTranscriptData, IsofoxTranscriptDataComparer, TestIsofoxTranscriptDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new IsofoxTranscriptDataComparer(new ComparConfig());
        builder = TestIsofoxTranscriptDataBuilder.BUILDER;
        IsofoxTranscriptData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                FLD_GENE_NAME, b -> b.geneName = alternateValueSource.TranscriptExpression().geneName(),
                FLD_ADJ_TPM, b -> b.tpm = alternateValueSource.TranscriptExpression().tpm()
        );
        nameToAlternateIndexInitializer = Map.of(
                FLD_TRANS_NAME, b -> b.transcriptName = alternateValueSource.TranscriptExpression().transcriptName()
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
    }

    @Override
    @Test
    public void fullyMatchesSelfInReportableMode()
    {
        // Override since Isofox output is never compared in reportable mode
    }

    @Override
    @Test
    public void singleFieldMismatchesAreRecognizedInReportableMode()
    {
        // Override since Isofox output is never compared in reportable mode
    }
}
