package com.hartwig.hmftools.isofox.novel.cohort;

import java.util.List;

import com.hartwig.hmftools.isofox.novel.AltSpliceJunction;

public class AltSjFilter
{
    public final List<String> RestrictedGeneIds;
    public final List<String> ExcludedGeneIds;
    public final int MinFragments;

    public AltSjFilter(final List<String> restrictedGeneIds, final List<String> excludedGeneIds, final int minFragments)
    {
        RestrictedGeneIds = restrictedGeneIds;
        ExcludedGeneIds = excludedGeneIds;
        MinFragments = minFragments;
    }

    public boolean passesFilter(final AltSpliceJunction altSJ)
    {
        if(MinFragments > 0 && altSJ.getFragmentCount() < MinFragments)
            return false;

        if(!RestrictedGeneIds.isEmpty() && !RestrictedGeneIds.contains(altSJ.getGeneId()))
            return false;

        if(!ExcludedGeneIds.isEmpty() && ExcludedGeneIds.contains(altSJ.getGeneId()))
            return false;

        return true;
    }
}
