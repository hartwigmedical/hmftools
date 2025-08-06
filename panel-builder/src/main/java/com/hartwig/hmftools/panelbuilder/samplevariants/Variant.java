package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAG_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_FRAG_COUNT_MIN_STRICT;

import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.wisp.CategoryType;

import org.jetbrains.annotations.Nullable;

public abstract class Variant
{
    private SelectionStatus mStatus;

    public Variant()
    {
        mStatus = SelectionStatus.NOT_SET;
    }

    public abstract CategoryType categoryType();

    public abstract String description();

    @Nullable
    public String gene()
    {
        return null;
    }

    public abstract double copyNumber();

    public abstract double vaf();

    @Nullable
    public String otherData()
    {
        return null;
    }

    public abstract int tumorFragments();

    public boolean hasPhaseVariants()
    {
        return false;
    }

    public abstract boolean reported();

    public abstract VariantProbeData generateProbe(final RefGenomeInterface refGenome);

    public abstract boolean checkFilters();

    public abstract boolean checkAndRegisterLocation(ProximateLocations registeredLocations);

    public boolean checkAndRegisterGeneLocation(Map<String, Integer> geneDisruptions)
    {
        return true;
    }

    public boolean passNonReportableFilters(boolean strictLimits)
    {
        return true;
    }

    protected static boolean passesFragmentCountLimit(int fragmentCount, boolean strictLimits)
    {
        int limit = strictLimits ? SAMPLE_FRAG_COUNT_MIN_STRICT : SAMPLE_FRAG_COUNT_MIN;
        return fragmentCount >= limit;
    }

    public SelectionStatus selectionStatus()
    {
        return mStatus;
    }

    public boolean isSelected()
    {
        return mStatus == SelectionStatus.SELECTED;
    }

    public void setSelectionStatus(final SelectionStatus status)
    {
        mStatus = status;
    }
}
