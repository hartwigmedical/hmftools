package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.purple.GeneCopyNumberData.FLD_MAX_COPY_NUMBER;
import static com.hartwig.hmftools.compar.purple.GeneCopyNumberData.FLD_MIN_COPY_NUMBER;

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
        comparer = new GeneCopyNumberComparer(new ComparConfig());
        builder = TestGeneCopyNumberDataBuilder.BUILDER;
        GeneCopyNumberData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                FLD_MIN_COPY_NUMBER, b -> b.minCopyNumber = alternateValueSource.CopyNumber.minCopyNumber(),
                FLD_MAX_COPY_NUMBER, b -> b.maxCopyNumber = alternateValueSource.CopyNumber.maxCopyNumber()
        );
        nameToAlternateIndexInitializer = Map.of("Gene", b -> b.gene = alternateValueSource.CopyNumber.geneName());
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
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
