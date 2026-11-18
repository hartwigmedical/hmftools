package com.hartwig.hmftools.compar.cider;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class Cdr3LocusSummaryDataTest extends ComparableItemTest<Cdr3LocusSummaryData, Cdr3LocusSummaryComparer, TestCdr3LocusSummaryDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new Cdr3LocusSummaryComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestCdr3LocusSummaryDataBuilder.BUILDER;
        Cdr3LocusSummaryData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer =
                Map.of(Cdr3LocusSummaryComparer.Fields.PassSequences.toString(), b -> b.passSequenceCount = alternateValueSource.Cdr3LocusSummary.passSequences());
        nameToAlternateIndexInitializer = Map.of("Locus", b -> b.locus = alternateValueSource.Cdr3LocusSummary.locus());
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
