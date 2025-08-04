package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_GC_MAX;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_GC_MAX_LOWER;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_GC_MIN;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_GC_MIN_LOWER;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.FRAG_COUNT_MIN;
import static com.hartwig.hmftools.panelbuilder.samplevariants.Constants.FRAG_COUNT_MIN_LOWER;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.wisp.CategoryType;

public abstract class Variant
{
    private String mSequence;
    private SelectionStatus mStatus;

    public Variant()
    {
        mSequence = "";
        mStatus = SelectionStatus.NOT_SET;
    }

    public abstract CategoryType categoryType();

    public abstract String description();

    public abstract String gene();

    public List<String> refSequences()
    {
        return Lists.newArrayList();
    }

    public void setSequence(final String sequence)
    {
        mSequence = sequence;
    }

    public String sequence()
    {
        return mSequence;
    }

    public abstract double copyNumber();

    public abstract double vaf();

    public double gc()
    {
        return calcGcPercent(sequence());
    }

    public String otherData()
    {
        return "";
    }

    public abstract int tumorFragments();

    public boolean hasPhaseVariants()
    {
        return false;
    }

    public abstract boolean reported();

    public abstract void generateSequences(final RefGenomeInterface refGenome);

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
        if(useLowerLimits)
        {
            return gcRatio >= SAMPLE_GC_MIN_LOWER && gcRatio <= SAMPLE_GC_MAX_LOWER;
        }
        else
        {
            return gcRatio >= SAMPLE_GC_MIN && gcRatio <= SAMPLE_GC_MAX;
        }
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

    public int sequenceCount()
    {
        return 1 + refSequences().size();
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
