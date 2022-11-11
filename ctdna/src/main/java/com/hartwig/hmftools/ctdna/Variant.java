package com.hartwig.hmftools.ctdna;

import static com.hartwig.hmftools.ctdna.SelectionStatus.NOT_SET;
import static com.hartwig.hmftools.ctdna.SelectionStatus.SELECTED;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public abstract class Variant
{
    private String mSequence;
    private SelectionStatus mStatus;

    public Variant()
    {
        mSequence = "";
        mStatus = NOT_SET;
    }

    abstract CategoryType categoryType();

    abstract String description();

    abstract String gene();

    public List<String> refSequences() { return Lists.newArrayList(); }

    public void setSequence(final String sequence) { mSequence = sequence; }
    public String sequence() { return mSequence; }

    abstract double copyNumber();

    abstract double vaf();

    public double gc() { return VariantUtils.calcGcPercent(sequence()); }

    String otherData() { return ""; }

    abstract int tumorFragments();

    boolean hasPhaseVariants() { return false; }

    abstract boolean reported();

    abstract void generateSequences(final RefGenomeInterface refGenome, final PvConfig config);

    abstract boolean checkAndRegisterLocation(final ProximateLocations registeredLocations);

    boolean passNonReportableFilters(final PvConfig config) { return true; }

    int sequenceCount() { return 1 + refSequences().size(); }

    public SelectionStatus selectionStatus() { return mStatus; }
    public boolean isSelected() { return mStatus == SELECTED; }
    public void setSelectionStatus(final SelectionStatus status) { mStatus = status; }
}
