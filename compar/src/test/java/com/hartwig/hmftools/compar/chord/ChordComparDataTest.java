package com.hartwig.hmftools.compar.chord;

import java.util.Collections;
import java.util.Map;

import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItemTest;

import org.junit.Before;

public class ChordComparDataTest extends ComparableItemTest<ChordComparData, ChordComparer, TestChordComparDataBuilder>
{
    @Before
    public void setUp()
    {
        comparer = new ChordComparer(new ComparConfig(), Collections.emptyMap());
        builder = TestChordComparDataBuilder.BUILDER;
        ChordComparData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                ChordComparer.Fields.BRCA1.toString(), b -> b.brca1 = alternateValueSource.Chord.BRCA1Value(),
                ChordComparer.Fields.BRCA2.toString(), b -> b.brca2 = alternateValueSource.Chord.BRCA2Value(),
                ChordComparer.Fields.Status.toString(), b -> b.status = alternateValueSource.Chord.hrStatus(),
                ChordComparer.Fields.Type.toString(), b -> b.type = alternateValueSource.Chord.hrdType(),
                ChordComparer.Fields.Score.toString(), b -> b.score = alternateValueSource.Chord.hrdValue()
        );
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
