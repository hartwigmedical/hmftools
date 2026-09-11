package com.hartwig.hmftools.compar.cider;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class CiderVdjDataTest extends ComparableItemTest<CiderVdjData, CiderVdjComparer, TestCiderVdjDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new CiderVdjComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestCiderVdjDataBuilder.BUILDER;
        CiderVdjData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                CiderVdjComparer.Fields.Filter.toString(), b -> b.filter = alternateValueSource.mCdr3Sequence.filter(),
                CiderVdjComparer.Fields.Locus.toString(), b -> b.locus = alternateValueSource.mCdr3Sequence.locus()
        );
        nameToAlternateIndexInitializer = Map.of("cdr3Seq", b -> b.cdr3Seq = alternateValueSource.mCdr3Sequence.cdr3Seq());
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
