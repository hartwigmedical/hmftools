package com.hartwig.hmftools.common.variant.structural;

import org.immutables.value.Value;

@Value.Immutable
public abstract class StructuralVariantData {

    public abstract int id();
    public abstract String startChromosome();
    public abstract String endChromosome();
    public abstract int startPosition();
    public abstract int endPosition();
    public abstract byte startOrientation();
    public abstract byte endOrientation();
    public abstract String startHomologySequence();
    public abstract String endHomologySequence();
    public abstract double startAF();
    public abstract double endAF();
    public abstract double ploidy();
    public abstract double adjustedStartAF();
    public abstract double adjustedEndAF();
    public abstract double adjustedStartCopyNumber();
    public abstract double adjustedEndCopyNumber();
    public abstract double adjustedStartCopyNumberChange();
    public abstract double adjustedEndCopyNumberChange();
    public abstract String insertSequence();
    public abstract StructuralVariantType type();
    public abstract String filter();
    public abstract boolean imprecise();
    public abstract double qualityScore();
    public abstract String event();
    public abstract int startTumorVariantFragmentCount();
    public abstract int startTumorReferenceFragmentCount();
    public abstract int startNormalVariantFragmentCount();
    public abstract int startNormalReferenceFragmentCount();
    public abstract int endTumorVariantFragmentCount();
    public abstract int endTumorReferenceFragmentCount();
    public abstract int endNormalVariantFragmentCount();
    public abstract int endNormalReferenceFragmentCount();
    public abstract int startIntervalOffsetStart();
    public abstract int startIntervalOffsetEnd();
    public abstract int endIntervalOffsetStart();
    public abstract int endIntervalOffsetEnd();
    public abstract int inexactHomologyOffsetStart();
    public abstract int inexactHomologyOffsetEnd();
    public abstract String startLinkedBy();
    public abstract String endLinkedBy();
    public abstract String vcfId();
    public abstract boolean recovered();
    public abstract String recoveryMethod();
    public abstract String recoveryFilter();
    public abstract String startRefContext();
    public abstract String endRefContext();
    public abstract String insertSequenceAlignments();
    public abstract String insertSequenceRepeatClass();
    public abstract String insertSequenceRepeatType();
    public abstract Byte insertSequenceRepeatOrientation();
    public abstract double insertSequenceRepeatCoverage();
    public abstract int startAnchoringSupportDistance();
    public abstract int endAnchoringSupportDistance();

}
