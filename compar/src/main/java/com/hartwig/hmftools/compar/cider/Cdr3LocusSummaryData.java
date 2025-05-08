package com.hartwig.hmftools.compar.cider;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.Category.CDR3_LOCUS_SUMMARY;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;

import static org.apache.commons.lang3.StringUtils.capitalize;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.cider.Cdr3LocusSummary;
import com.hartwig.hmftools.common.cider.Cdr3LocusSummaryFile;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class Cdr3LocusSummaryData implements ComparableItem
{
    public final Cdr3LocusSummary Cdr3LocusSummary;

    public Cdr3LocusSummaryData(final Cdr3LocusSummary cdr3LocusSummary)
    {
        Cdr3LocusSummary = cdr3LocusSummary;
    }

    @Override
    public Category category() { return CDR3_LOCUS_SUMMARY; }

    @Override
    public String key()
    {
        return String.format("locus(%s)", Cdr3LocusSummary.locus());
    }

    public static List<String> comparedFieldNames()
    {
        return List.of(capitalize(Cdr3LocusSummaryFile.Column.passSequences.name()));
    }

    @Override
    public List<String> displayValues()
    {
        return List.of(String.valueOf(Cdr3LocusSummary.passSequences()));
    }

    @Override
    public boolean reportable() { return true; }

    @Override
    public boolean isPass() {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem comparableItem)
    {
        final Cdr3LocusSummary other = ((Cdr3LocusSummaryData) comparableItem).Cdr3LocusSummary;
        return Cdr3LocusSummary.locus().equals(other.locus());
    }

    @Override
    public Mismatch findMismatch(final ComparableItem comparableItem, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        final Cdr3LocusSummary other = ((Cdr3LocusSummaryData) comparableItem).Cdr3LocusSummary;

        final List<String> diffs = new ArrayList<>();
        checkDiff(diffs, Cdr3LocusSummaryFile.Column.passSequences.name(), Cdr3LocusSummary.passSequences(), other.passSequences(), thresholds);

        return createMismatchFromDiffs(this, comparableItem, diffs, matchLevel, includeMatches);
    }

    public String toString()
    {
        return format("locus(%s) pass sequences(%d)", Cdr3LocusSummary.locus(), Cdr3LocusSummary.passSequences());
    }
}
