package com.hartwig.hmftools.isofox.novel.cohort;

import java.util.List;

import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;

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

    public boolean passesFilter(final AltSpliceJunctionFile altSJ)
    {
        if(MinFragments > 0 && altSJ.FragmentCount < MinFragments)
            return false;

        if(!RestrictedGeneIds.isEmpty() && !RestrictedGeneIds.contains(altSJ.GeneId))
            return false;

        if(!ExcludedGeneIds.isEmpty() && ExcludedGeneIds.contains(altSJ.GeneId))
            return false;

        return true;
    }
}
