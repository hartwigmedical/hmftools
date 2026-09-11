package com.hartwig.hmftools.compar.cider;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.CDR3_SEQUENCE;

import java.util.List;

import com.hartwig.hmftools.common.cider.Cdr3Sequence;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class CiderVdjData extends ComparableItem
{
    public final Cdr3Sequence mCdr3Sequence;

    public CiderVdjData(final Cdr3Sequence cdr3Sequence, final List<FieldInfo> fields)
    {
        mCdr3Sequence = cdr3Sequence;

        addStringValue(CiderVdjComparer.Fields.Filter.toString(), cdr3Sequence.filter().toString(), fields);
        addStringValue(CiderVdjComparer.Fields.Locus.toString(), cdr3Sequence.locus().toString(), fields);
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
