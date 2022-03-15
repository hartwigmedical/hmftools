package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.classification.CkbEventAndGeneExtractor;
import com.hartwig.hmftools.ckb.classification.CkbProteinAnnotationExtractor;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.immuno.ActionableHLA;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.extraction.ActionableEventFactory;
import com.hartwig.hmftools.serve.extraction.EventExtractor;
import com.hartwig.hmftools.serve.extraction.EventExtractorOutput;
import com.hartwig.hmftools.serve.extraction.ExtractionFunctions;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.ImmutableExtractionResult;
import com.hartwig.hmftools.serve.extraction.codon.CodonAnnotation;
import com.hartwig.hmftools.serve.extraction.codon.CodonFunctions;
import com.hartwig.hmftools.serve.extraction.codon.ImmutableKnownCodon;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberFunctions;
import com.hartwig.hmftools.serve.extraction.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.events.EventInterpretation;
import com.hartwig.hmftools.serve.extraction.events.ImmutableEventInterpretation;
import com.hartwig.hmftools.serve.extraction.exon.ExonAnnotation;
import com.hartwig.hmftools.serve.extraction.exon.ExonFunctions;
import com.hartwig.hmftools.serve.extraction.exon.ImmutableKnownExon;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.FusionFunctions;
import com.hartwig.hmftools.serve.extraction.fusion.ImmutableKnownFusionPair;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.extraction.hotspot.ImmutableKnownHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.util.ProgressTracker;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CkbExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CkbExtractor.class);

    @NotNull
    private final EventExtractor eventExtractor;

    public CkbExtractor(@NotNull final EventExtractor eventExtractor) {
        this.eventExtractor = eventExtractor;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<CkbEntry> ckbEntries) {
        List<ExtractionResult> extractions = Lists.newArrayList();
        List<EventInterpretation> interpretation = Lists.newArrayList();

        ProgressTracker tracker = new ProgressTracker("CKB", ckbEntries.size());
        for (CkbEntry entry : ckbEntries) {
            // Assume entries without variants are filtered out prior to extraction
            assert !entry.variants().isEmpty();

            Variant variant = entry.variants().get(0);
            String gene = CkbEventAndGeneExtractor.extractGene(variant);
            String event = CkbEventAndGeneExtractor.extractEvent(variant);

            if (entry.type() == EventType.UNKNOWN) {
                LOGGER.warn("No event type known for '{}' on '{}'", event, gene);
            }

            String transcript = null;

            EventExtractorOutput eventExtractorOutput = eventExtractor.extract(gene, transcript, entry.type(), event, null);
            Set<? extends ActionableEvent> actionableEvents = ActionableEntryFactory.toActionableEntries(entry, variant.variant());

            interpretation.add(ImmutableEventInterpretation.builder()
                    .source(Knowledgebase.CKB)
                    .sourceEvent(variant.variant())
                    .interpretGene(gene)
                    .interpretEvent(event)
                    .interpretEventType(entry.type())
                    .build());
            extractions.add(toExtractionResult(gene, event, transcript, eventExtractorOutput, actionableEvents, interpretation));

            tracker.update();
        }

        return ExtractionFunctions.merge(extractions);
    }

    @NotNull
    private static ExtractionResult toExtractionResult(@NotNull String gene, @NotNull String variant, @Nullable String transcript,
            @NotNull EventExtractorOutput output, @NotNull Set<? extends ActionableEvent> actionableEvents,
            @NotNull List<EventInterpretation> interpretation) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();
        Set<ActionableHLA> actionableHLA = Sets.newHashSet();

        for (ActionableEvent event : actionableEvents) {
            actionableHotspots.addAll(ActionableEventFactory.toActionableHotspots(event, output.hotspots()));
            actionableRanges.addAll(ActionableEventFactory.toActionableRanges(event, output.codons()));
            actionableRanges.addAll(ActionableEventFactory.toActionableRanges(event, output.exons()));

            if (output.geneLevelEvent() != null) {
                actionableGenes.add(ActionableEventFactory.geneLevelEventToActionableGene(event, output.geneLevelEvent()));
            }

            if (output.knownCopyNumber() != null) {
                actionableGenes.add(ActionableEventFactory.copyNumberToActionableGene(event, output.knownCopyNumber()));
            }

            if (output.knownFusionPair() != null) {
                actionableFusions.add(ActionableEventFactory.toActionableFusion(event, output.knownFusionPair()));
            }

            if (output.characteristic() != null) {
                actionableCharacteristics.add(ActionableEventFactory.toActionableCharacteristic(event, output.characteristic()));
            }

            if (output.hla() != null) {
                actionableHLA.add(ActionableEventFactory.toActionableHLa(event, output.hla()));
            }
        }

        return ImmutableExtractionResult.builder()
                .eventInterpretation(interpretation)
                .refGenomeVersion(Knowledgebase.CKB.refGenomeVersion())
                .knownHotspots(convertToKnownHotspots(output.hotspots(), gene, variant, transcript))
                .knownCodons(convertToKnownCodons(output.codons()))
                .knownExons(convertToKnownExons(output.exons()))
                .knownCopyNumbers(convertToKnownAmpsDels(output.knownCopyNumber()))
                .knownFusionPairs(convertToKnownFusions(output.knownFusionPair()))
                .actionableHotspots(actionableHotspots)
                .actionableRanges(actionableRanges)
                .actionableGenes(actionableGenes)
                .actionableFusions(actionableFusions)
                .actionableCharacteristics(actionableCharacteristics)
                .actionableHLA(actionableHLA)
                .build();
    }

    @NotNull
    private static Set<KnownHotspot> convertToKnownHotspots(@Nullable List<VariantHotspot> hotspots, @NotNull String gene,
            @NotNull String variant, @Nullable String transcript) {
        Set<KnownHotspot> knownHotspots = Sets.newHashSet();

        if (hotspots != null) {
            CkbProteinAnnotationExtractor proteinExtractor = new CkbProteinAnnotationExtractor();
            for (VariantHotspot hotspot : hotspots) {
                knownHotspots.add(ImmutableKnownHotspot.builder()
                        .from(hotspot)
                        .addSources(Knowledgebase.CKB)
                        .gene(gene)
                        .transcript(transcript)
                        .proteinAnnotation(proteinExtractor.apply(variant))
                        .build());
            }
        }

        return HotspotFunctions.consolidate(knownHotspots);
    }

    @NotNull
    private static Set<KnownCodon> convertToKnownCodons(@Nullable List<CodonAnnotation> codonAnnotations) {
        Set<KnownCodon> codons = Sets.newHashSet();

        if (codonAnnotations != null) {
            for (CodonAnnotation codonAnnotation : codonAnnotations) {
                codons.add(ImmutableKnownCodon.builder().annotation(codonAnnotation).addSources(Knowledgebase.CKB).build());
            }
        }
        return CodonFunctions.consolidate(codons);
    }

    @NotNull
    private static Set<KnownExon> convertToKnownExons(@Nullable List<ExonAnnotation> exonAnnotations) {
        Set<KnownExon> exons = Sets.newHashSet();

        if (exonAnnotations != null) {
            for (ExonAnnotation exonAnnotation : exonAnnotations) {
                exons.add(ImmutableKnownExon.builder().annotation(exonAnnotation).addSources(Knowledgebase.CKB).build());
            }
        }
        return ExonFunctions.consolidate(exons);
    }

    @NotNull
    private static Set<KnownCopyNumber> convertToKnownAmpsDels(@Nullable KnownCopyNumber knownCopyNumber) {
        Set<KnownCopyNumber> copyNumbers = Sets.newHashSet();
        if (knownCopyNumber != null) {
            copyNumbers.add(ImmutableKnownCopyNumber.builder().from(knownCopyNumber).addSources(Knowledgebase.CKB).build());
        }
        return CopyNumberFunctions.consolidate(copyNumbers);
    }

    @NotNull
    private static Set<KnownFusionPair> convertToKnownFusions(@Nullable KnownFusionPair knownFusionPair) {
        Set<KnownFusionPair> fusions = Sets.newHashSet();
        if (knownFusionPair != null) {
            fusions.add(ImmutableKnownFusionPair.builder().from(knownFusionPair).addSources(Knowledgebase.CKB).build());
        }
        return FusionFunctions.consolidate(fusions);
    }
}