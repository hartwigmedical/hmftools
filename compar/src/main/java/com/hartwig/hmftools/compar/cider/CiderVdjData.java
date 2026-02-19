package com.hartwig.hmftools.compar.cider;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.CDR3_SEQUENCE;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import static org.apache.commons.lang3.StringUtils.capitalize;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cider.Cdr3Sequence;
import com.hartwig.hmftools.common.cider.Cdr3SequenceFile;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class CiderVdjData implements ComparableItem
{
    public final Cdr3Sequence mCdr3Sequence;

    public static final String FILTER_FIELD = capitalize(Cdr3SequenceFile.Column.filter.name());
    public static final String LOCUS_FIELD = capitalize(Cdr3SequenceFile.Column.locus.name());

    public CiderVdjData(final Cdr3Sequence mCdr3Sequence)
    {
        this.mCdr3Sequence = mCdr3Sequence;
    }

    @Override
    public CategoryType category() { return CDR3_SEQUENCE; }

    @Override
    public String key()
    {
        return String.format("cdr3AA(%s) cdr3Seq(%s)", mCdr3Sequence.cdr3AA(), mCdr3Sequence.cdr3Seq());
    }

    public static List<String> comparedFieldNames()
    {
        return List.of(FILTER_FIELD, LOCUS_FIELD);
    }

    @Override
    public List<String> displayValues()
    {
        return List.of(mCdr3Sequence.filter(), mCdr3Sequence.locus());
    }

    @Override
    public boolean matches(final ComparableItem o)
    {
        final Cdr3Sequence other = ((CiderVdjData) o).mCdr3Sequence;

        // match cdr3 seq by their sequence
        return mCdr3Sequence.cdr3Seq().equals(other.cdr3Seq());
    }

    @Override
    public Mismatch findMismatch(final ComparableItem o, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final Cdr3Sequence other = ((CiderVdjData) o).mCdr3Sequence;

        final List<String> diffs = new ArrayList<>();

        checkDiff(diffs, FILTER_FIELD, mCdr3Sequence.filter(), other.filter());
        checkDiff(diffs, LOCUS_FIELD, mCdr3Sequence.locus(), other.locus());

        return createMismatchFromDiffs(this, o, diffs, matchLevel, includeMatches);
    }

    public String toString()
    {
        return format("cdr3AA(%s)", mCdr3Sequence.cdr3AA());
    }
}
