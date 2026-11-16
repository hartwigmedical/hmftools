package com.hartwig.hmftools.compar.isofox;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class RnaSummaryData extends ComparableItem
{
    public final RnaStatistics Statistics;

    public RnaSummaryData(final RnaStatistics statistics, final List<FieldInfo> fields)
    {
        Statistics = statistics;

        addStringValue(
                RnaSummaryComparer.Fields.QcStatus.toString(), Statistics.qcStatus().stream()
                        .map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM)),
                fields);

        addLongValue(RnaSummaryComparer.Fields.TotalFragments.toString(), Statistics.totalFragments(), fields);
        addLongValue(RnaSummaryComparer.Fields.DuplicateFragments.toString(), Statistics.duplicateFragments(), fields);

        addDoubleValue(RnaSummaryComparer.Fields.SplicedFragmentPerc.toString(), Statistics.splicedFragmentPerc(), fields);
        addDoubleValue(RnaSummaryComparer.Fields.UnsplicedFragmentPerc.toString(), Statistics.unsplicedFragmentPerc(), fields);
        addDoubleValue(RnaSummaryComparer.Fields.AltFragmentPerc.toString(), Statistics.altFragmentPerc(), fields);
        addDoubleValue(RnaSummaryComparer.Fields.ChimericFragmentPerc.toString(), Statistics.chimericFragmentPerc(), fields);
        addIntValue(RnaSummaryComparer.Fields.SplicedGeneCount.toString(), Statistics.splicedGeneCount(), fields);
        addDoubleValue(RnaSummaryComparer.Fields.FragLength5th.toString(), Statistics.fragmentLength5thPercent(), fields);
        addDoubleValue(RnaSummaryComparer.Fields.FragLength50th.toString(), Statistics.fragmentLength50thPercent(), fields);
        addDoubleValue(RnaSummaryComparer.Fields.FragLength95th.toString(), Statistics.fragmentLength95thPercent(), fields);
        addDoubleValue(RnaSummaryComparer.Fields.MedianGCRatio.toString(), Statistics.medianGCRatio(), fields);
        addDoubleValue(RnaSummaryComparer.Fields.ForwardStrandPercent.toString(), Statistics.forwardStrandPercent(), fields);
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
