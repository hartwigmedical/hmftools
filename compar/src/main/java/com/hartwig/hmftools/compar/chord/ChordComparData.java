package com.hartwig.hmftools.compar.chord;

import static com.hartwig.hmftools.compar.common.CategoryType.CHORD;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class ChordComparData implements ComparableItem
{
    public final ChordData Chord;

    protected static final String FLD_BRCA1 = "BRCA1";
    protected static final String FLD_BRCA2 = "BRCA2";
    protected static final String FLD_STATUS = "Status";
    protected static final String FLD_TYPE = "Type";
    protected static final String FLD_SCORE = "Score";

    public ChordComparData(final ChordData chord)
    {
        Chord = chord;
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
