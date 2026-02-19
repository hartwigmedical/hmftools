package com.hartwig.hmftools.compar.chord;

import static com.hartwig.hmftools.compar.chord.ChordComparData.FLD_BRCA1;
import static com.hartwig.hmftools.compar.chord.ChordComparData.FLD_BRCA2;
import static com.hartwig.hmftools.compar.chord.ChordComparData.FLD_SCORE;
import static com.hartwig.hmftools.compar.chord.ChordComparData.FLD_STATUS;
import static com.hartwig.hmftools.compar.chord.ChordComparData.FLD_TYPE;

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
        comparer = new ChordComparer(new ComparConfig());
        builder = TestChordComparDataBuilder.BUILDER;
        ChordComparData alternateValueSource = builder.createWithAlternateDefaults();
        fieldToAlternateValueInitializer = Map.of(
                FLD_BRCA1, b -> b.brca1 = alternateValueSource.Chord.BRCA1Value(),
                FLD_BRCA2, b -> b.brca2 = alternateValueSource.Chord.BRCA2Value(),
                FLD_STATUS, b -> b.status = alternateValueSource.Chord.hrStatus(),
                FLD_TYPE, b -> b.type = alternateValueSource.Chord.hrdType(),
                FLD_SCORE, b -> b.score = alternateValueSource.Chord.hrdValue()
        );
        nameToAlternateIndexInitializer = Collections.emptyMap();
        reportabilityFieldToFalseReportabilityInitializer = Collections.emptyMap();
        nameToNonPassInitializer = Collections.emptyMap();
    }
}
