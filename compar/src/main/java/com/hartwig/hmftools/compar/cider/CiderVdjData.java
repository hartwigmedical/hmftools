package com.hartwig.hmftools.compar.cider;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.CDR3_SEQUENCE;

import static org.apache.commons.lang3.StringUtils.capitalize;

import com.hartwig.hmftools.common.cider.Cdr3Sequence;
import com.hartwig.hmftools.common.cider.Cdr3SequenceFile;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

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

    @Override
    public boolean matches(final ComparableItem o)
    {
        final Cdr3Sequence other = ((CiderVdjData) o).mCdr3Sequence;

        // match cdr3 seq by their sequence
        return mCdr3Sequence.cdr3Seq().equals(other.cdr3Seq());
    }

    public String toString()
    {
        return format("cdr3AA(%s)", mCdr3Sequence.cdr3AA());
    }
}
