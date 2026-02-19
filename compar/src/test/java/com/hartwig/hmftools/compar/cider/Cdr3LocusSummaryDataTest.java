package com.hartwig.hmftools.compar.cider;

import static com.hartwig.hmftools.compar.cider.Cdr3LocusSummaryData.PASS_SEQUENCES_FIELD;

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
        comparer = new Cdr3LocusSummaryComparer(new ComparConfig());
        builder = TestCdr3LocusSummaryDataBuilder.BUILDER;
        Cdr3LocusSummaryData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer =
                Map.of(PASS_SEQUENCES_FIELD, b -> b.passSequenceCount = alternateValueSource.Cdr3LocusSummary.passSequences());
        nameToAlternateIndexInitializer = Map.of("Locus", b -> b.locus = alternateValueSource.Cdr3LocusSummary.locus());
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
