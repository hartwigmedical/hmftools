package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.classification.CkbConstants;
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
import com.hartwig.hmftools.serve.extraction.codon.ImmutableCodonAnnotation;
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
    private static final String DELIMITER = ",";
    @NotNull
    private final EventExtractor eventExtractor;

    public CkbExtractor(@NotNull final EventExtractor eventExtractor) {
        this.eventExtractor = eventExtractor;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<CkbEntry> ckbEntries) {
        List<ExtractionResult> extractions = Lists.newArrayList();

        ProgressTracker tracker = new ProgressTracker("CKB", ckbEntries.size());
        for (CkbEntry entry : ckbEntries) {
            // Assume entries without variants are filtered out prior to extraction
            assert !entry.variants().isEmpty();

            int variantCount = entry.variants().size();
            Variant variant = entry.variants().get(0);
            String event = variantCount > 1 ? concat(entry.variants()) : CkbEventAndGeneExtractor.extractEvent(variant);
            String gene = variantCount > 1 ? "Multiple" : CkbEventAndGeneExtractor.extractGene(variant);

            if (entry.type() == EventType.UNKNOWN) {
                LOGGER.warn("No event type known for '{}' on '{}'", event, gene);
            } else {
                EventExtractorOutput extraction = eventExtractor.extract(gene, null, entry.type(), event);
                String sourceEvent;
                if (!gene.equals(CkbConstants.NO_GENE)) {
                    sourceEvent = gene + " " + event;
                } else {
                    sourceEvent = event;
                }

                Set<ActionableEntry> actionableEvents = ActionableEntryFactory.toActionableEntries(entry, sourceEvent);

                EventInterpretation interpretation = ImmutableEventInterpretation.builder()
                        .source(Knowledgebase.CKB)
                        .sourceEvent(sourceEvent)
                        .interpretedGene(gene)
                        .interpretedEvent(event)
                        .interpretedEventType(entry.type())
                        .build();

                extractions.add(toExtractionResult(gene, event, null, extraction, actionableEvents, interpretation));
            }

            tracker.update();
        }

        return ExtractionFunctions.merge(extractions);
    }

    @NotNull
    private static String concat(@NotNull List<Variant> variants) {
        StringJoiner joiner = new StringJoiner(DELIMITER);
        for (Variant variant : variants) {
            joiner.add(variant.variant());
        }
        return joiner.toString();
    }

    @NotNull
    private static ExtractionResult toExtractionResult(@NotNull String gene, @NotNull String variant, @Nullable String transcript,
            @NotNull EventExtractorOutput output, @NotNull Set<ActionableEntry> actionableEvents,
            @NotNull EventInterpretation interpretation) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();
        Set<ActionableHLA> actionableHLA = Sets.newHashSet();

        List<CodonAnnotation> codons = Lists.newArrayList();

        for (ActionableEvent event : actionableEvents) {
            codons = curateCodons(output.codons());

            actionableHotspots.addAll(ActionableEventFactory.toActionableHotspots(event, output.hotspots()));
            actionableRanges.addAll(ActionableEventFactory.toActionableRanges(event, codons));
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
                .refGenomeVersion(Knowledgebase.CKB.refGenomeVersion())
                .eventInterpretations(Lists.newArrayList(interpretation))
                .knownHotspots(convertToKnownHotspots(output.hotspots(), gene, variant, transcript))
                .knownCodons(convertToKnownCodons(codons))
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

    @VisibleForTesting
    @Nullable
    static List<CodonAnnotation> curateCodons(@Nullable List<CodonAnnotation> codonAnnotations) {
        if (codonAnnotations == null) {
            return null;
        }

        List<CodonAnnotation> curatedCodons = Lists.newArrayList();
        for (CodonAnnotation codon : codonAnnotations) {
            if (codon.gene().equals("BRAF") && codon.rank() == 600) {
                curatedCodons.add(ImmutableCodonAnnotation.builder()
                        .from(codon)
                        .transcript("ENST00000646891")
                        .start(140753335)
                        .end(140753337)
                        .build());
            } else {
                curatedCodons.add(codon);
            }
        }
        return curatedCodons;
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
