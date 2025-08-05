package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.FRAG_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.FRAG_COUNT_MIN_STRICT;

import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.wisp.CategoryType;
import com.hartwig.hmftools.panelbuilder.PanelCoverage;
import com.hartwig.hmftools.panelbuilder.ProbeFactory;
import com.hartwig.hmftools.panelbuilder.ProbeGenerationResult;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;

import org.jetbrains.annotations.Nullable;

public abstract class Variant
{
    // TODO: refactor this, also figure out GC constraints
    @Nullable
    private ProbeGenerationResult mProbeGenResult;
    private SelectionStatus mStatus;

    public Variant()
    {
        mProbeGenResult = null;
        mStatus = SelectionStatus.NOT_SET;
    }

    public abstract CategoryType categoryType();

    public abstract String description();

    @Nullable
    public String gene()
    {
        return null;
    }

    protected void setProbeGenResult(final ProbeGenerationResult result)
    {
        mProbeGenResult = result;
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

    public abstract void generateProbe(final RefGenomeInterface refGenome, final ProbeFactory probeFactory, final PanelCoverage coverage);

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
        int limit = strictLimits ? FRAG_COUNT_MIN_STRICT : FRAG_COUNT_MIN;
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

    protected TargetMetadata targetMetadata()
    {
        return new TargetMetadata(TargetMetadata.Type.SAMPLE_VARIANT, description());
    }
}
