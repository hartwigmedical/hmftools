package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.rna.RnaStatisticFile.Column;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.rna.RnaQcFilter;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

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
        return CategoryType.RNA_SUMMARY;
    }

    @Override
    public String key()
    {
        return "";
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

    public static String qcStatus(final List<RnaQcFilter> status)
    {
        StringJoiner sj = new StringJoiner(";");
        status.forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }
}
