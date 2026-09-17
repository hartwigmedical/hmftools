package com.hartwig.hmftools.compar.sigs;

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
        comparer = new SigsComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestSigsDataBuilder.BUILDER;
        SigsData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                SigsComparer.Fields.Percent.toString(), b -> b.percent = alternateValueSource.Allocation.percent()
        );
        nameToAlternateIndexInitializer = Map.of(
                "Signature", b -> b.signature = alternateValueSource.Allocation.signature()
        );
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }

    @Override
    @Test
    public void fullyMatchesSelfInReportableMode()
    {
        // Override since signatures output is never compared in reportable mode
    }

    @Override
    @Test
    public void singleFieldMismatchesAreRecognizedInReportableMode()
    {
        // Override since signatures output is never compared in reportable mode
    }
}