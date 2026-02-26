package com.hartwig.hmftools.compar.sigs;

import static com.hartwig.hmftools.common.sigs.SignatureAllocationFile.PERCENT_FLD;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;
import org.junit.Test;

public class SigsDataTest extends ComparableItemTest<SigsData, SigsComparer, TestSigsDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new SigsComparer(new ComparConfig());
        builder = TestSigsDataBuilder.BUILDER;
        SigsData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                PERCENT_FLD, b -> b.percent = alternateValueSource.SignatureAllocation().percent()
        );
        nameToAlternateIndexInitializer = Map.of(
                "Signature", b -> b.signature = alternateValueSource.SignatureAllocation().signature()
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