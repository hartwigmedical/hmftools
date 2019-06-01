package com.hartwig.hmftools.common.variant.structural;

import org.immutables.value.Value;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
public abstract class StructuralVariantData {

    @Nullable
    public abstract String id();
    public abstract String startChromosome();
    public abstract String endChromosome();
    public abstract long startPosition();
    public abstract long endPosition();
    public abstract byte startOrientation();
    public abstract byte endOrientation();
    public abstract String startHomologySequence();
    public abstract String endHomologySequence();
    public abstract Double startAF();
    public abstract Double endAF();
    public abstract Double ploidy();
    public abstract Double adjustedStartAF();
    public abstract Double adjustedEndAF();
    public abstract Double adjustedStartCopyNumber();
    public abstract Double adjustedEndCopyNumber();
    public abstract Double adjustedStartCopyNumberChange();
    public abstract Double adjustedEndCopyNumberChange();
    public abstract String insertSequence();
    public abstract StructuralVariantType type();
    public abstract String filter();
    public abstract boolean imprecise();
    public abstract Double qualityScore();
    public abstract String event();
    public abstract Integer startTumorVariantFragmentCount();
    public abstract Integer startTumorReferenceFragmentCount();
    public abstract Integer startNormalVariantFragmentCount();
    public abstract Integer startNormalReferenceFragmentCount();
    public abstract Integer endTumorVariantFragmentCount();
    public abstract Integer endTumorReferenceFragmentCount();
    public abstract Integer endNormalVariantFragmentCount();
    public abstract Integer endNormalReferenceFragmentCount();
    public abstract Integer startIntervalOffsetStart();
    public abstract Integer startIntervalOffsetEnd();
    public abstract Integer endIntervalOffsetStart();
    public abstract Integer endIntervalOffsetEnd();
    public abstract Integer inexactHomologyOffsetStart();
    public abstract Integer inexactHomologyOffsetEnd();
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
    public abstract Double insertSequenceRepeatCoverage();
    public abstract int startAnchoringSupportDistance();
    public abstract int endAnchoringSupportDistance();

}
