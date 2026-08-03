package com.hartwig.hmftools.compar.purple;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;
import org.junit.Test;

public class GeneCopyNumberDataTest extends ComparableItemTest<GeneCopyNumberData, GeneCopyNumberComparer, TestGeneCopyNumberDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new GeneCopyNumberComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestGeneCopyNumberDataBuilder.BUILDER;
        GeneCopyNumberData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                GeneCopyNumberComparer.Fields.MinCopyNumber.toString(), b -> b.minCopyNumber = alternateValueSource.CopyNumber.minCopyNumber(),
                GeneCopyNumberComparer.Fields.MaxCopyNumber.toString(), b -> b.maxCopyNumber = alternateValueSource.CopyNumber.maxCopyNumber()
        );
        nameToAlternateIndexInitializer = Map.of("Gene", b -> b.gene = alternateValueSource.CopyNumber.geneName());
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }

    @Override
    @Test
    public void fullyMatchesSelfInReportableMode()
    {
        // Override since copy numbers are never compared in reportable mode
    }

    @Override
    @Test
    public void singleFieldMismatchesAreRecognizedInReportableMode()
    {
        // Override since copy numbers are never compared in reportable mode
    }
}
