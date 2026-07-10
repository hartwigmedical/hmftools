package com.hartwig.hmftools.compar.isofox;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.rna.RnaQcFilter;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.CategoryType;

public record IsofoxSummaryData(RnaStatistics RnaStatistics) implements ComparableItem
{
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
