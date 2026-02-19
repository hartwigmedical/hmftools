package com.hartwig.hmftools.compar.cider;

import static com.hartwig.hmftools.compar.cider.CiderVdjData.FILTER_FIELD;
import static com.hartwig.hmftools.compar.cider.CiderVdjData.LOCUS_FIELD;

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
        comparer = new CiderVdjComparer(new ComparConfig());
        builder = TestCiderVdjDataBuilder.BUILDER;
        CiderVdjData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                FILTER_FIELD, b -> b.filter = alternateValueSource.mCdr3Sequence.filter(),
                LOCUS_FIELD, b -> b.locus = alternateValueSource.mCdr3Sequence.locus()
        );
        nameToAlternateIndexInitializer = Map.of("cdr3Seq", b -> b.cdr3Seq = alternateValueSource.mCdr3Sequence.cdr3Seq());
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
