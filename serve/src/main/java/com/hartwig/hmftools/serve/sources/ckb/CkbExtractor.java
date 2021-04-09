package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.common.refseq.RefSeq;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
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
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CkbExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CkbExtractor.class);

    @NotNull
    private final EventExtractor eventExtractor;
    @NotNull
    private final ActionableEvidenceFactory actionableEvidenceFactory;

    public CkbExtractor(@NotNull final EventExtractor eventExtractor, @NotNull final ActionableEvidenceFactory actionableEvidenceFactory) {
        this.eventExtractor = eventExtractor;
        this.actionableEvidenceFactory = actionableEvidenceFactory;
    }

    @Nullable
    public String extractCanonicalTranscript(@NotNull String refseqToMatch, @NotNull List<RefSeq> refSeqMatchFile) {
        for (RefSeq refSeq: refSeqMatchFile) {
            if (refSeq.dbPrimaryAcc().equals(refseqToMatch)) {
                return refSeq.transcriptId();
            }
        }
        return null;
    }

    @NotNull
    public ExtractionResult extract(@NotNull List<CkbEntry> ckbEntries, @NotNull List<RefSeq> refSeqMatchFile) {
        List<ExtractionResult> extractions = Lists.newArrayList();
        List<EventExtractorOutput> eventExtractions = Lists.newArrayList();

        ProgressTracker tracker = new ProgressTracker("CKB", ckbEntries.size());
        for (CkbEntry entry : ckbEntries) {
            if (entry.variants().size() == 1) {
                Variant variant = entry.variants().get(0);
                String event;
                String geneSymbol;

                if (variant.variant().equals("fusion") && variant.impact() != null && variant.impact().equals("fusion")) {
                    event = "fusion promiscuous";
                    geneSymbol = variant.gene().geneSymbol();
                } else if (variant.impact() != null && variant.impact().equals("fusion")) {
                    event = variant.variant().replaceAll("\\s+", "");
                    geneSymbol = variant.fullName();
                } else if (variant.variant().contains("exon")) {
                    event = variant.variant().replace("exon", "exon ");
                    geneSymbol = variant.gene().geneSymbol();
                } else {
                    event = variant.variant();
                    geneSymbol = variant.gene().geneSymbol();
                }

                eventExtractions.add(eventExtractor.extract(geneSymbol,
                        extractCanonicalTranscript(entry.variants().get(0).gene().canonicalTranscript(), refSeqMatchFile),
                        entry.type(),
                        event));
                Set<ActionableEvent> actionableEvents = actionableEvidenceFactory.toActionableEvents(entry);

                extractions.add(toExtractionResult(actionableEvents, eventExtractions));
                extractions.add(ImmutableExtractionResult.builder()
                        .knownHotspots(convertToHotspots(eventExtractions, entry))
                        .knownCodons(convertToCodons(eventExtractions))
                        .knownExons(convertToExons(eventExtractions))
                        .knownCopyNumbers(convertToKnownAmpsDels(eventExtractions))
                        .knownFusionPairs(convertToKnownFusions(eventExtractions))
                        .build());

                if (entry.type() == EventType.UNKNOWN) {
                    LOGGER.warn("No event type known for '{}' on '{}'",
                            entry.variants().get(0).variant(),
                            entry.variants().get(0).gene().geneSymbol());
                }
            }
            tracker.update();
        }

        actionableEvidenceFactory.evaluateCuration();

        return ExtractionFunctions.merge(extractions);
    }

    @NotNull
    private static Set<KnownHotspot> convertToHotspots(@NotNull List<EventExtractorOutput> results, @NotNull CkbEntry entry) {
        ProteinAnnotationExtractor proteinExtractor = new ProteinAnnotationExtractor();
        Set<KnownHotspot> hotspots = Sets.newHashSet();
        for (EventExtractorOutput result : results) {
            if (result.hotspots() != null) {
                for (VariantHotspot hotspot : result.hotspots()) {
                    hotspots.add(ImmutableKnownHotspot.builder()
                            .from(hotspot)
                            .addSources(Knowledgebase.CKB)
                            .gene(entry.variants().get(0).gene().geneSymbol())
                            .transcript(entry.variants().get(0).gene().canonicalTranscript())
                            .proteinAnnotation(proteinExtractor.apply(entry.variants().get(0).variant()))
                            .build());
                }
            }
        }
        return HotspotFunctions.consolidate(hotspots);
    }

    @NotNull
    private static Set<KnownCodon> convertToCodons(@NotNull List<EventExtractorOutput> results) {
        Set<KnownCodon> codons = Sets.newHashSet();

        for (EventExtractorOutput result : results) {
            if (result.codons() != null) {
                for (CodonAnnotation codonAnnotation : result.codons()) {
                    codons.add(ImmutableKnownCodon.builder().annotation(codonAnnotation).addSources(Knowledgebase.CKB).build());
                }
            }
        }
        return CodonFunctions.consolidate(codons);
    }

    @NotNull
    private static Set<KnownExon> convertToExons(@NotNull List<EventExtractorOutput> results) {
        Set<KnownExon> exons = Sets.newHashSet();

        for (EventExtractorOutput result : results) {
            if (result.exons() != null) {
                for (ExonAnnotation exonAnnotation : result.exons()) {
                    exons.add(ImmutableKnownExon.builder().annotation(exonAnnotation).addSources(Knowledgebase.CKB).build());
                }
            }
        }
        return ExonFunctions.consolidate(exons);
    }

    @NotNull
    private static Set<KnownCopyNumber> convertToKnownAmpsDels(@NotNull List<EventExtractorOutput> results) {
        Set<KnownCopyNumber> copyNumbers = Sets.newHashSet();
        for (EventExtractorOutput result : results) {
            if (result.knownCopyNumber() != null) {
                copyNumbers.add(ImmutableKnownCopyNumber.builder().from(result.knownCopyNumber()).addSources(Knowledgebase.CKB).build());

            }
        }
        return CopyNumberFunctions.consolidate(copyNumbers);
    }

    @NotNull
    private static Set<KnownFusionPair> convertToKnownFusions(@NotNull List<EventExtractorOutput> results) {
        Set<KnownFusionPair> fusions = Sets.newHashSet();
        for (EventExtractorOutput result : results) {
            if (result.knownFusionPair() != null) {
                fusions.add(ImmutableKnownFusionPair.builder().from(result.knownFusionPair()).addSources(Knowledgebase.CKB).build());
            }
        }
        return FusionFunctions.consolidate(fusions);
    }

    @NotNull
    private static ExtractionResult toExtractionResult(@NotNull Set<ActionableEvent> actionableTrials,
            @NotNull List<EventExtractorOutput> eventExtractions) {
        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        Set<ActionableGene> actionableGenes = Sets.newHashSet();
        Set<ActionableFusion> actionableFusions = Sets.newHashSet();
        Set<ActionableCharacteristic> actionableCharacteristics = Sets.newHashSet();

        for (ActionableEvent trial : actionableTrials) {
            for (EventExtractorOutput extraction : eventExtractions) {
                actionableHotspots.addAll(ActionableEventFactory.toActionableHotspots(trial, extraction.hotspots()));
                actionableRanges.addAll(ActionableEventFactory.toActionableRanges(trial, extraction.codons()));
                actionableRanges.addAll(ActionableEventFactory.toActionableRanges(trial, extraction.exons()));

                if (extraction.geneLevelEvent() != null) {
                    actionableGenes.add(ActionableEventFactory.geneLevelEventToActionableGene(trial, extraction.geneLevelEvent()));
                }

                if (extraction.knownCopyNumber() != null) {
                    actionableGenes.add(ActionableEventFactory.copyNumberToActionableGene(trial, extraction.knownCopyNumber()));
                }

                if (extraction.knownFusionPair() != null) {
                    actionableFusions.add(ActionableEventFactory.toActionableFusion(trial, extraction.knownFusionPair()));
                }

                if (extraction.characteristic() != null) {
                    actionableCharacteristics.add(ActionableEventFactory.toActionableCharacteristic(trial, extraction.characteristic()));
                }
            }
        }

        return ImmutableExtractionResult.builder()
                .actionableHotspots(actionableHotspots)
                .actionableRanges(actionableRanges)
                .actionableGenes(actionableGenes)
                .actionableFusions(actionableFusions)
                .actionableCharacteristics(actionableCharacteristics)
                .build();
    }
}