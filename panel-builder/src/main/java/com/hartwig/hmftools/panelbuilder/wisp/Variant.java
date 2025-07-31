package com.hartwig.hmftools.panelbuilder.wisp;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;

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

    abstract CategoryType categoryType();

    abstract String description();

    abstract String gene();

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

    abstract double copyNumber();

    abstract double vaf();

    public double gc()
    {
        return calcGcPercent(sequence());
    }

    String otherData()
    {
        return "";
    }

    abstract int tumorFragments();

    boolean hasPhaseVariants()
    {
        return false;
    }

    abstract boolean reported();

    abstract void generateSequences(final RefGenomeInterface refGenome, final ProbeConfig config);

    abstract boolean checkFilters();

    abstract boolean checkAndRegisterLocation(final ProximateLocations registeredLocations);

    boolean checkAndRegisterGeneLocation(final Map<String, Integer> geneDisruptions)
    {
        return true;
    }

    boolean passNonReportableFilters(final ProbeConfig config, boolean useLowerLimits)
    {
        return true;
    }

    protected static boolean passesGcRatioLimit(double gcRatio, final ProbeConfig config, boolean useLowerLimits)
    {
        if(useLowerLimits)
        {
            return gcRatio >= config.GcRatioLimitLowerMin && gcRatio <= config.GcRatioLimitLowerMax;
        }
        else
        {
            return gcRatio >= config.GcRatioLimitMin && gcRatio <= config.GcRatioLimitMax;
        }
    }

    protected static boolean passesFragmentCountLimit(int fragmentCount, final ProbeConfig config, boolean useLowerLimits)
    {
        if(useLowerLimits)
        {
            return fragmentCount >= config.FragmentCountMinLower;
        }
        else
        {
            return fragmentCount >= config.FragmentCountMin;
        }
    }

    int sequenceCount()
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
