package com.hartwig.hmftools.compar.purple;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;
import org.junit.Test;

public class CopyNumberDataTest extends ComparableItemTest<CopyNumberData, CopyNumberComparer, TestCopyNumberDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new CopyNumberComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestCopyNumberDataBuilder.BUILDER;
        CopyNumberData alternateValueSource = builder.createWithAlternateDefaults();

        fieldToAlternateValueInitializer = Map.of(
                CopyNumberComparer.Fields.CopyNumber.toString(), b -> b.copyNumber = alternateValueSource.CopyNumber,
                CopyNumberComparer.Fields.MajorAlleleCopyNumber.toString(), b -> b.majorAlleleCopyNumber = alternateValueSource.MajorAlleleCopyNumber,
                CopyNumberComparer.Fields.Method.toString(), b -> b.method = alternateValueSource.Method
        );

        nameToAlternateIndexInitializer = Map.of("Chromosome", b -> b.chromosome = alternateValueSource.Chromosome);
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
