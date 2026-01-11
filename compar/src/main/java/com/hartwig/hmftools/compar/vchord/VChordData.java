package com.hartwig.hmftools.compar.vchord;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.V_CHORD;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.vchord.VChordPrediction;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public record VChordData(VChordPrediction VChord) implements ComparableItem
{
    static final String FLD_BREAST = "BreastCancerHrdScore";
    static final String FLD_OVARIAN = "OvarianCancerHrdScore";
    static final String FLD_PANCREATIC = "PancreaticCancerScore";
    static final String FLD_PROSTATE = "ProstateCancerScore";
    static final String FLD_OTHER = "OtherCancerScore";

    @Override
    public CategoryType category()
    {
        return V_CHORD;
    }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%.2f", VChord.breastCancerHrdScore()));
        values.add(format("%.2f", VChord.ovarianCancerHrdScore()));
        values.add(format("%.2f", VChord.pancreaticCancerScore()));
        values.add(format("%.2f", VChord.prostateCancerScore()));
        values.add(format("%.2f", VChord.otherCancerScore()));
        return values;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final VChordData otherData = (VChordData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_BREAST, VChord.breastCancerHrdScore(), otherData.VChord.breastCancerHrdScore(), thresholds);
        checkDiff(diffs, FLD_OVARIAN, VChord.ovarianCancerHrdScore(), otherData.VChord.ovarianCancerHrdScore(), thresholds);
        checkDiff(diffs, FLD_PANCREATIC, VChord.pancreaticCancerScore(), otherData.VChord.pancreaticCancerScore(), thresholds);
        checkDiff(diffs, FLD_PROSTATE, VChord.prostateCancerScore(), otherData.VChord.prostateCancerScore(), thresholds);
        checkDiff(diffs, FLD_OTHER, VChord.otherCancerScore(), otherData.VChord.otherCancerScore(), thresholds);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }
}
