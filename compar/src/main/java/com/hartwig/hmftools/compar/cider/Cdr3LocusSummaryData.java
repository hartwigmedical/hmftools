package com.hartwig.hmftools.compar.cider;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.CategoryType.CDR3_LOCUS_SUMMARY;

import com.hartwig.hmftools.common.cider.Cdr3LocusSummary;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public class Cdr3LocusSummaryData implements ComparableItem
{
    public final Cdr3LocusSummary Cdr3LocusSummary;

    public Cdr3LocusSummaryData(final Cdr3LocusSummary cdr3LocusSummary)
    {
        Cdr3LocusSummary = cdr3LocusSummary;
    }

    @Override
    public CategoryType category() { return CDR3_LOCUS_SUMMARY; }

    @Override
    public String key()
    {
        return String.format("locus(%s)", Cdr3LocusSummary.locus());
    }

    @Override
    public boolean matches(final ComparableItem comparableItem)
    {
        final Cdr3LocusSummary other = ((Cdr3LocusSummaryData) comparableItem).Cdr3LocusSummary;
        return Cdr3LocusSummary.locus().equals(other.locus());
    }

    public String toString()
    {
        return format("locus(%s) pass sequences(%d)", Cdr3LocusSummary.locus(), Cdr3LocusSummary.passSequences());
    }
}
