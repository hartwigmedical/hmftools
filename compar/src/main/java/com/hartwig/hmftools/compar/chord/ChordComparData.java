package com.hartwig.hmftools.compar.chord;

import static com.hartwig.hmftools.compar.common.CategoryType.CHORD;

import java.util.List;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class ChordComparData extends ComparableItem
{
    public final ChordData Chord;

    public ChordComparData(final ChordData chord, final List<FieldInfo> fields)
    {
        Chord = chord;

        addDoubleValue(ChordComparer.Fields.BRCA1.toString(), chord.BRCA1Value(), fields);
        addDoubleValue(ChordComparer.Fields.BRCA2.toString(), chord.BRCA2Value(), fields);
        addDoubleValue(ChordComparer.Fields.Score.toString(), chord.hrdValue(), fields);
        addStringValue(ChordComparer.Fields.Type.toString(), chord.hrdType(), fields);
        addStringValue(ChordComparer.Fields.Status.toString(), chord.hrStatus().toString(), fields);
    }

    @Override
    public CategoryType category() { return CHORD; }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }
}
