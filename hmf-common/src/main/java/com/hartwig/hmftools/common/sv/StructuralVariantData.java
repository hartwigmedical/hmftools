package com.hartwig.hmftools.common.sv;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

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
    public abstract double junctionCopyNumber();
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

    public static StructuralVariantData convertSvData(final EnrichedStructuralVariant var, int svId)
    {
        return ImmutableStructuralVariantData.builder()
                .id(svId)
                .startChromosome(var.chromosome(true))
                .endChromosome(var.end() == null ? "0" : var.chromosome(false))
                .startPosition(var.position(true).intValue())
                .endPosition(var.end() == null ? -1 : var.position(false).intValue())
                .startOrientation(var.orientation(true))
                .endOrientation(var.end() == null ? (byte) 0 : var.orientation(false))
                .startHomologySequence(var.start().homology())
                .endHomologySequence(var.end() == null ? "" : var.end().homology())
                .junctionCopyNumber(getValueNotNull(var.junctionCopyNumber()))
                .startAF(getValueNotNull(var.start().alleleFrequency()))
                .endAF(var.end() == null ? 0 : getValueNotNull(var.end().alleleFrequency()))
                .adjustedStartAF(getValueNotNull(var.start().adjustedAlleleFrequency()))
                .adjustedEndAF(var.end() == null ? 0 : getValueNotNull(var.end().adjustedAlleleFrequency()))
                .adjustedStartCopyNumber(getValueNotNull(var.start().adjustedCopyNumber()))
                .adjustedEndCopyNumber(var.end() == null ? 0 : getValueNotNull(var.end().adjustedCopyNumber()))
                .adjustedStartCopyNumberChange(getValueNotNull(var.start().adjustedCopyNumberChange()))
                .adjustedEndCopyNumberChange(var.end() == null ? 0 : getValueNotNull(var.end().adjustedCopyNumberChange()))
                .insertSequence(var.insertSequence())
                .type(var.type())
                .filter(var.filter())
                .imprecise(var.imprecise())
                .qualityScore(getValueNotNull(var.qualityScore()))
                .event(getValueNotNull(var.event()))
                .startTumorVariantFragmentCount(getValueNotNull(var.start().tumorVariantFragmentCount()))
                .startTumorReferenceFragmentCount(getValueNotNull(var.start().tumorReferenceFragmentCount()))
                .startNormalVariantFragmentCount(getValueNotNull(var.start().normalVariantFragmentCount()))
                .startNormalReferenceFragmentCount(getValueNotNull(var.start().normalReferenceFragmentCount()))
                .endTumorVariantFragmentCount(var.end() == null ? 0 : getValueNotNull(var.end().tumorVariantFragmentCount()))
                .endTumorReferenceFragmentCount(var.end() == null ? 0 : getValueNotNull(var.end().tumorReferenceFragmentCount()))
                .endNormalVariantFragmentCount(var.end() == null ? 0 : getValueNotNull(var.end().normalVariantFragmentCount()))
                .endNormalReferenceFragmentCount(var.end() == null ? 0 : getValueNotNull(var.end().normalReferenceFragmentCount()))
                .startIntervalOffsetStart(getValueNotNull(var.start().startOffset()))
                .startIntervalOffsetEnd(getValueNotNull(var.start().endOffset()))
                .endIntervalOffsetStart(var.end() == null ? 0 : getValueNotNull(var.end().startOffset()))
                .endIntervalOffsetEnd(var.end() == null ? 0 : getValueNotNull(var.end().endOffset()))
                .inexactHomologyOffsetStart(getValueNotNull(var.start().inexactHomologyOffsetStart()))
                .inexactHomologyOffsetEnd(getValueNotNull(var.start().inexactHomologyOffsetEnd()))
                .startLinkedBy(getValueNotNull(var.startLinkedBy()))
                .endLinkedBy(getValueNotNull(var.endLinkedBy()))
                .vcfId(getValueNotNull(var.id()))
                .startRefContext(getValueNotNull(var.start().refGenomeContext()))
                .endRefContext(var.end() == null ? "" : getValueNotNull(var.end().refGenomeContext()))
                .recovered(var.recovered())
                .recoveryMethod((getValueNotNull(var.recoveryMethod())))
                .recoveryFilter(getValueNotNull(var.recoveryFilter()))
                .insertSequenceAlignments(getValueNotNull(var.insertSequenceAlignments()))
                .insertSequenceRepeatClass(getValueNotNull(var.insertSequenceRepeatClass()))
                .insertSequenceRepeatType(getValueNotNull(var.insertSequenceRepeatType()))
                .insertSequenceRepeatOrientation(getValueNotNull(var.insertSequenceRepeatOrientation()))
                .insertSequenceRepeatCoverage(getValueNotNull(var.insertSequenceRepeatCoverage()))
                .startAnchoringSupportDistance(var.start().anchoringSupportDistance())
                .endAnchoringSupportDistance(var.end() == null ? 0 : var.end().anchoringSupportDistance())
                .build();
    }

    public static double getValueNotNull(@Nullable Double value) {
        return value != null ? value : 0D;
    }

    public static int getValueNotNull(@Nullable Integer value) {
        return value != null ? value : 0;
    }

    public static byte getValueNotNull(@Nullable Byte value) {
        return value != null ? value : 0;
    }

    @NotNull
    public static String getValueNotNull(@Nullable String value) {
        return value != null ? value : Strings.EMPTY;
    }


}
