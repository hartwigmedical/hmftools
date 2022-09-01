package com.hartwig.hmftools.serve;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.serve.actionability.ImmutableTreatment;
import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.characteristic.ImmutableActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.immuno.ActionableHLA;
import com.hartwig.hmftools.serve.actionability.immuno.ImmutableActionableHLA;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.actionability.range.RangeType;
import com.hartwig.hmftools.serve.cancertype.ImmutableCancerType;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableCodonAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableKnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.extraction.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableExonAnnotation;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableKnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.extraction.hotspot.ImmutableKnownHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.extraction.util.MutationTypeFilter;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ServeTestFactory {

    private ServeTestFactory() {
    }

    @NotNull
    public static ExtractionResult createResultForSource(@NotNull Knowledgebase source) {
        return ImmutableExtractionResult.builder()
                .refGenomeVersion(source.refGenomeVersion())
                .addKnownHotspots(createTestKnownHotspotForSource(source))
                .addKnownCodons(createTestKnownCodonForSource(source))
                .addKnownExons(createTestKnownExonForSource(source))
                .addKnownCopyNumbers(createTestKnownCopyNumberForSource(source))
                .addKnownFusionPairs(createTestKnownFusionPairForSource(source))
                .addActionableHotspots(createTestActionableHotspotForSource(source))
                .addActionableRanges(createTestActionableRangeForSource(source))
                .addActionableGenes(createTestActionableGeneForSource(source))
                .addActionableFusions(createTestActionableFusionForSource(source))
                .addActionableCharacteristics(createTestActionableCharacteristicForSource(source))
                .addActionableHLA(createTestActionableImmunoHLAForSource(source))
                .build();
    }

    @NotNull
    public static KnownHotspot createTestKnownHotspotForSource(@NotNull Knowledgebase source) {
        return ImmutableKnownHotspot.builder().from(createTestKnownHotspot()).addSources(source).build();
    }

    @NotNull
    public static KnownHotspot createTestKnownHotspot() {
        return ImmutableKnownHotspot.builder()
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .proteinAnnotation(Strings.EMPTY)
                .build();
    }

    @NotNull
    public static KnownCodon createTestKnownCodonForSource(@NotNull Knowledgebase source) {
        return ImmutableKnownCodon.builder().from(createTestKnownCodon()).addSources(source).build();
    }

    @NotNull
    public static KnownCodon createTestKnownCodon() {
        return ImmutableKnownCodon.builder()
                .annotation(ImmutableCodonAnnotation.builder().from(createTestCodonAnnotation()).build())
                .build();
    }

    @NotNull
    public static CodonAnnotation createTestCodonAnnotation() {
        return ImmutableCodonAnnotation.builder()
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .start(0)
                .end(0)
                .mutationType(MutationTypeFilter.ANY)
                .rank(0)
                .build();
    }

    @NotNull
    public static KnownExon createTestKnownExonForSource(@NotNull Knowledgebase source) {
        return ImmutableKnownExon.builder().from(createTestKnownExon()).addSources(source).build();
    }

    @NotNull
    public static KnownExon createTestKnownExon() {
        return ImmutableKnownExon.builder().annotation(ImmutableExonAnnotation.builder().from(createTestExonAnnotation()).build()).build();
    }

    @NotNull
    public static ExonAnnotation createTestExonAnnotation() {
        return ImmutableExonAnnotation.builder()
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .start(0)
                .end(0)
                .mutationType(MutationTypeFilter.ANY)
                .rank(0)
                .build();
    }

    @NotNull
    public static KnownCopyNumber createTestKnownCopyNumberForSource(@NotNull Knowledgebase source) {
        return ImmutableKnownCopyNumber.builder().from(createTestKnownCopyNumber()).addSources(source).build();
    }

    @NotNull
    public static KnownCopyNumber createTestKnownCopyNumber() {
        return ImmutableKnownCopyNumber.builder().gene(Strings.EMPTY).type(CopyNumberType.AMPLIFICATION).build();
    }

    @NotNull
    public static KnownFusionPair createTestKnownFusionPairForSource(@NotNull Knowledgebase source) {
        return ImmutableKnownFusionPair.builder().from(createTestKnownFusionPair()).addSources(source).build();
    }

    @NotNull
    public static KnownFusionPair createTestKnownFusionPair() {
        return ImmutableKnownFusionPair.builder().geneUp(Strings.EMPTY).geneDown(Strings.EMPTY).build();
    }

    @NotNull
    public static ActionableHotspot createTestActionableHotspotForSource(@NotNull Knowledgebase source) {
        return ImmutableActionableHotspot.builder().from(createTestActionableHotspot()).source(source).build();
    }

    @NotNull
    public static ActionableHotspot createTestActionableHotspot() {
        return ImmutableActionableHotspot.builder()
                .from(createTestBaseEvent())
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .build();
    }

    @NotNull
    public static ActionableRange createTestActionableRangeForSource(@NotNull Knowledgebase source) {
        return ImmutableActionableRange.builder().from(createTestActionableRange()).source(source).build();
    }

    @NotNull
    public static ActionableRange createTestActionableRange() {
        return ImmutableActionableRange.builder()
                .from(createTestBaseEvent())
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .start(0)
                .end(0)
                .mutationType(MutationTypeFilter.ANY)
                .rangeType(RangeType.EXON)
                .rank(0)
                .build();
    }

    @NotNull
    public static ActionableGene createTestActionableGeneForSource(@NotNull Knowledgebase source) {
        return ImmutableActionableGene.builder().from(createTestActionableGene()).source(source).build();
    }

    @NotNull
    public static ActionableGene createTestActionableGene() {
        return ImmutableActionableGene.builder()
                .from(createTestBaseEvent())
                .gene(Strings.EMPTY)
                .event(GeneLevelEvent.AMPLIFICATION)
                .build();
    }

    @NotNull
    public static ActionableFusion createTestActionableFusionForSource(@NotNull Knowledgebase source) {
        return ImmutableActionableFusion.builder().from(createTestActionableFusion()).source(source).build();
    }

    @NotNull
    public static ActionableFusion createTestActionableFusion() {
        return ImmutableActionableFusion.builder().from(createTestBaseEvent()).geneUp(Strings.EMPTY).geneDown(Strings.EMPTY).build();
    }

    @NotNull
    public static ActionableCharacteristic createTestActionableCharacteristicForSource(@NotNull Knowledgebase source) {
        return ImmutableActionableCharacteristic.builder().from(createTestActionableCharacteristic()).source(source).build();
    }

    @NotNull
    public static ActionableCharacteristic createTestActionableCharacteristic() {
        return ImmutableActionableCharacteristic.builder()
                .from(createTestBaseEvent())
                .name(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE)
                .build();
    }

    @NotNull
    public static ActionableHLA createTestActionableImmunoHLAForSource(@NotNull Knowledgebase source) {
        return ImmutableActionableHLA.builder().from(createTestActionableHLA()).source(source).build();
    }

    @NotNull
    public static ActionableHLA createTestActionableHLA() {
        return ImmutableActionableHLA.builder().from(createTestBaseEvent()).hlaType(Strings.EMPTY).build();
    }

    @NotNull
    private static ActionableEvent createTestBaseEvent() {
        return createTestBaseEvent(Knowledgebase.HARTWIG_CURATED);
    }

    @NotNull
    private static ActionableEvent createTestBaseEvent(@NotNull Knowledgebase source) {
        return ActionabilityTestUtil.create(source,
                "source event",
                Sets.newHashSet(),
                ImmutableTreatment.builder()
                        .treament("treatment")
                        .sourceRelevantTreatmentApproaches(Sets.newHashSet("drugClasses"))
                        .relevantTreatmentApproaches(Sets.newHashSet("drugClasses"))
                        .build(),
                ImmutableCancerType.builder().name("applicable name").doid("applicable doid").build(),
                Sets.newHashSet(ImmutableCancerType.builder().name("blacklist name").doid("blacklist doid").build()),
                EvidenceLevel.A,
                EvidenceDirection.RESPONSIVE,
                Sets.newHashSet());
    }
}