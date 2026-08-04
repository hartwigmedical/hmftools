package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_TRANS_NAME;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;
import org.junit.Test;

public class RnaTranscriptDataTest extends ComparableItemTest<RnaTranscriptData, RnaTranscriptDataComparer, TestIsofoxTranscriptDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new RnaTranscriptDataComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestIsofoxTranscriptDataBuilder.BUILDER;
        RnaTranscriptData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                RnaTranscriptDataComparer.Fields.AdjTPM.toString(), b -> b.tpm = alternateValueSource.Expression.tpm()
        );
        nameToAlternateIndexInitializer = Map.of(
                FLD_TRANS_NAME, b -> b.transcriptName = alternateValueSource.Expression.transcriptName()
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
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
