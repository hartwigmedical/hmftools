package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class IsofoxSummaryData extends ComparableItem
{
    public final RnaStatistics Statistics;

    public IsofoxSummaryData(final RnaStatistics statistics, final List<FieldInfo> fields)
    {
        Statistics = statistics;

        addStringValue(
                IsofoxSummaryComparer.Fields.QcStatus.toString(), Statistics.qcStatus().stream()
                        .map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)),
                fields);

        addLongValue(IsofoxSummaryComparer.Fields.TotalFragments.toString(), Statistics.totalFragments(), fields);
        addLongValue(IsofoxSummaryComparer.Fields.DuplicateFragments.toString(), Statistics.duplicateFragments(), fields);

        addDoubleValue(IsofoxSummaryComparer.Fields.SplicedFragmentPerc.toString(), Statistics.splicedFragmentPerc(), fields);
        addDoubleValue(IsofoxSummaryComparer.Fields.UnsplicedFragmentPerc.toString(), Statistics.unsplicedFragmentPerc(), fields);
        addDoubleValue(IsofoxSummaryComparer.Fields.AltFragmentPerc.toString(), Statistics.altFragmentPerc(), fields);
        addDoubleValue(IsofoxSummaryComparer.Fields.ChimericFragmentPerc.toString(), Statistics.chimericFragmentPerc(), fields);
        addIntValue(IsofoxSummaryComparer.Fields.SplicedGeneCount.toString(), Statistics.splicedGeneCount(), fields);
        addDoubleValue(IsofoxSummaryComparer.Fields.FragLength5th.toString(), Statistics.fragmentLength5thPercent(), fields);
        addDoubleValue(IsofoxSummaryComparer.Fields.FragLength50th.toString(), Statistics.fragmentLength50thPercent(), fields);
        addDoubleValue(IsofoxSummaryComparer.Fields.FragLength95th.toString(), Statistics.fragmentLength95thPercent(), fields);
        addDoubleValue(IsofoxSummaryComparer.Fields.MedianGCRatio.toString(), Statistics.medianGCRatio(), fields);
        addDoubleValue(IsofoxSummaryComparer.Fields.ForwardStrandPercent.toString(), Statistics.forwardStrandPercent(), fields);
    }

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
}
