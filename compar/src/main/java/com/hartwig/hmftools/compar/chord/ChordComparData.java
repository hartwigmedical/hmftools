package com.hartwig.hmftools.compar.chord;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.Category.CHORD;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

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
    public Category category() { return CHORD; }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%.2f", Chord.BRCA1Value()));
        values.add(format("%.2f", Chord.BRCA2Value()));
        values.add(format("%.2f", Chord.hrdValue()));
        values.add(format("%s", Chord.hrStatus()));
        values.add(format("%s", Chord.hrdType()));
        return values;
    }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final ChordComparData otherData = (ChordComparData)other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_BRCA1, Chord.BRCA1Value(), otherData.Chord.BRCA1Value(), thresholds);
        checkDiff(diffs, FLD_BRCA2, Chord.BRCA2Value(), otherData.Chord.BRCA2Value(), thresholds);
        checkDiff(diffs, FLD_SCORE, Chord.hrdValue(), otherData.Chord.hrdValue(), thresholds);
        checkDiff(diffs, FLD_TYPE, Chord.hrdType(), otherData.Chord.hrdType());
        checkDiff(diffs, FLD_STATUS, Chord.hrStatus().toString(), otherData.Chord.hrStatus().toString());

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }
}
