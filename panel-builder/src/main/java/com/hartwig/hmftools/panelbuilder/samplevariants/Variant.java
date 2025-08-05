package com.hartwig.hmftools.panelbuilder.samplevariants;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_GC_TARGET;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_GC_TOLERANCE;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_GC_TOLERANCE_STRICT;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.FRAG_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.FRAG_COUNT_MIN_LOWER;

import java.util.Map;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.wisp.CategoryType;
import com.hartwig.hmftools.panelbuilder.PanelCoverage;
import com.hartwig.hmftools.panelbuilder.Probe;
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

    public Probe probe()
    {
        return mProbeGenResult.probes().get(0);
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

    public abstract boolean checkAndRegisterLocation(final ProximateLocations registeredLocations);

    public boolean checkAndRegisterGeneLocation(final Map<String, Integer> geneDisruptions)
    {
        return true;
    }

    public boolean passNonReportableFilters(boolean useLowerLimits)
    {
        return true;
    }

    protected static boolean passesGcRatioLimit(double gcRatio, boolean useLowerLimits)
    {
        double tolerance = useLowerLimits ? SAMPLE_GC_TOLERANCE : SAMPLE_GC_TOLERANCE_STRICT;
        return abs(gcRatio - SAMPLE_GC_TARGET) <= tolerance;
    }

    protected static boolean passesFragmentCountLimit(int fragmentCount, boolean useLowerLimits)
    {
        if(useLowerLimits)
        {
            return fragmentCount >= FRAG_COUNT_MIN_LOWER;
        }
        else
        {
            return fragmentCount >= FRAG_COUNT_MIN;
        }
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
