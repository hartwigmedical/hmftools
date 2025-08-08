package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAG_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAG_COUNT_MIN_STRICT;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.wisp.CategoryType;

public abstract class Variant
{
    public abstract CategoryType categoryType();

    public abstract double copyNumber();

    public abstract double vaf();

    public abstract int tumorFragments();

    public abstract boolean reported();

    public abstract VariantProbeData generateProbe(final RefGenomeInterface refGenome);

    public abstract boolean checkFilters();

    public abstract List<ProximateLocations.Location> checkedLocations();

    public boolean passNonReportableFilters(boolean strictLimits)
    {
        return true;
    }

    protected static boolean passesFragmentCountLimit(int fragmentCount, boolean strictLimits)
    {
        int limit = strictLimits ? SAMPLE_FRAG_COUNT_MIN_STRICT : SAMPLE_FRAG_COUNT_MIN;
        return fragmentCount >= limit;
    }
}
