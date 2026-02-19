package com.hartwig.hmftools.isofox.novel.cohort;

import java.util.List;

import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFile;

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

    public boolean passesFilter(final String geneId,  int fragCount)
    {
        if(MinFragments > 0 && fragCount < MinFragments)
            return false;

        if(!RestrictedGeneIds.isEmpty() && !RestrictedGeneIds.contains(geneId))
            return false;

        if(!ExcludedGeneIds.isEmpty() && ExcludedGeneIds.contains(geneId))
            return false;

        return true;
    }

    public boolean passesFilter(final AltSpliceJunctionFile altSJ)
    {
        return passesFilter(altSJ.GeneId, altSJ.FragmentCount);
    }
}
