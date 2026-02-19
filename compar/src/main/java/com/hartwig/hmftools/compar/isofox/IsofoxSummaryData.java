package com.hartwig.hmftools.compar.isofox;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.rna.RnaStatisticFile.Column;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.RnaQcFilter;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public record IsofoxSummaryData(RnaStatistics RnaStatistics) implements ComparableItem
{
    static final String FLD_QC_STATUS = Column.QcStatus.toString();
    static final String FLD_TOTAL_FRAGS = Column.TotalFragments.toString();
    static final String FLD_DUPLICATE_FRAGS = Column.DuplicateFragments.toString();
    static final String FLD_SPLICED_FRAG_PERC = Column.SplicedFragmentPerc.toString();
    static final String FLD_UNSPLICED_FRAG_PERC = Column.UnsplicedFragmentPerc.toString();
    static final String FLD_ALT_FRAG_PERC = Column.AltFragmentPerc.toString();
    static final String FLD_CHIMERIC_FRAG_PERC = Column.ChimericFragmentPerc.toString();
    static final String FLD_SPLICED_GENE_COUNT = Column.SplicedGeneCount.toString();
    static final String FLD_READ_LENGTH = Column.ReadLength.toString();
    static final String FLD_FRAG_LENGTH_5TH = Column.FragLength5th.toString();
    static final String FLD_FRAG_LENGTH_50TH = Column.FragLength50th.toString();
    static final String FLD_FRAG_LENGTH_95TH = Column.FragLength95th.toString();
    static final String FLD_ENRICHED_GENE_PERC = Column.EnrichedGenePercent.toString();
    static final String FLD_MEDIAN_GC_RATIO = Column.MedianGCRatio.toString();
    static final String FLD_FORWARD_STRAND_PERC = Column.ForwardStrandPercent.toString();

    @Override
    public CategoryType category()
    {
        return CategoryType.ISOFOX_SUMMARY;
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
        values.add(format("%s", qcStatus(RnaStatistics.qcStatus())));
        values.add(format("%d", RnaStatistics.totalFragments()));
        values.add(format("%d", RnaStatistics.duplicateFragments()));
        values.add(format("%.2f", RnaStatistics.splicedFragmentPerc()));
        values.add(format("%.2f", RnaStatistics.unsplicedFragmentPerc()));
        values.add(format("%.2f", RnaStatistics.altFragmentPerc()));
        values.add(format("%.2f", RnaStatistics.chimericFragmentPerc()));
        values.add(format("%d", RnaStatistics.readLength()));
        values.add(format("%.1f", RnaStatistics.fragmentLength5thPercent()));
        values.add(format("%.1f", RnaStatistics.fragmentLength50thPercent()));
        values.add(format("%.1f", RnaStatistics.fragmentLength95thPercent()));
        values.add(format("%.2f", RnaStatistics.enrichedGenePercent()));
        return values;
    }

    @Override
    public boolean reportable()
    {
        return false;
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
        final RnaStatistics ref = RnaStatistics;
        final RnaStatistics otherData = ((IsofoxSummaryData) other).RnaStatistics;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_QC_STATUS, qcStatus(ref.qcStatus()), qcStatus(otherData.qcStatus()));
        checkDiff(diffs, FLD_TOTAL_FRAGS, ref.totalFragments(), otherData.totalFragments(), thresholds);
        checkDiff(diffs, FLD_DUPLICATE_FRAGS, ref.duplicateFragments(), otherData.duplicateFragments(), thresholds);
        checkDiff(diffs, FLD_SPLICED_FRAG_PERC, ref.splicedFragmentPerc(), otherData.splicedFragmentPerc(), thresholds);
        checkDiff(diffs, FLD_UNSPLICED_FRAG_PERC, ref.unsplicedFragmentPerc(), otherData.unsplicedFragmentPerc(), thresholds);
        checkDiff(diffs, FLD_ALT_FRAG_PERC, ref.altFragmentPerc(), otherData.altFragmentPerc(), thresholds);
        checkDiff(diffs, FLD_CHIMERIC_FRAG_PERC, ref.chimericFragmentPerc(), otherData.chimericFragmentPerc(), thresholds);
        checkDiff(diffs, FLD_SPLICED_GENE_COUNT, ref.splicedGeneCount(), otherData.splicedGeneCount(), thresholds);
        checkDiff(diffs, FLD_READ_LENGTH, ref.readLength(), otherData.readLength());
        checkDiff(diffs, FLD_FRAG_LENGTH_5TH, ref.fragmentLength5thPercent(), otherData.fragmentLength5thPercent(), thresholds);
        checkDiff(diffs, FLD_FRAG_LENGTH_50TH, ref.fragmentLength50thPercent(), otherData.fragmentLength50thPercent(), thresholds);
        checkDiff(diffs, FLD_FRAG_LENGTH_95TH, ref.fragmentLength95thPercent(), otherData.fragmentLength95thPercent(), thresholds);
        checkDiff(diffs, FLD_ENRICHED_GENE_PERC, ref.enrichedGenePercent(), otherData.enrichedGenePercent(), thresholds);
        checkDiff(diffs, FLD_MEDIAN_GC_RATIO, ref.medianGCRatio(), otherData.medianGCRatio(), thresholds);
        checkDiff(diffs, FLD_FORWARD_STRAND_PERC, ref.forwardStrandPercent(), otherData.forwardStrandPercent(), thresholds);

        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }

    public static String qcStatus(final List<RnaQcFilter> status)
    {
        StringJoiner sj = new StringJoiner(";");
        status.forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }
}
