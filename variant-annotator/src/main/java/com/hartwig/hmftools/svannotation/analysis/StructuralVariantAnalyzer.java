package com.hartwig.hmftools.svannotation.analysis;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusionData;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusionModel;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svannotation.VariantAnnotator;
import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.GeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.GeneFusion;
import com.hartwig.hmftools.svannotation.annotations.ImmutableGeneDisruption;
import com.hartwig.hmftools.svannotation.annotations.ImmutableGeneFusion;
import com.hartwig.hmftools.svannotation.annotations.StructuralVariantAnnotation;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class StructuralVariantAnalyzer {

    @NotNull
    private final VariantAnnotator annotator;
    @NotNull
    private final Collection<HmfGenomeRegion> hmfGenePanelRegions;
    @NotNull
    private final CosmicFusionModel cosmicFusionModel;

    public StructuralVariantAnalyzer(@NotNull final VariantAnnotator annotator,
            @NotNull final Collection<HmfGenomeRegion> hmfGenePanelRegions, @NotNull final CosmicFusionModel cosmicFusionModel) {
        this.annotator = annotator;
        this.hmfGenePanelRegions = hmfGenePanelRegions;
        this.cosmicFusionModel = cosmicFusionModel;
    }

    @NotNull
    public StructuralVariantAnalysis run(@NotNull final List<StructuralVariant> variants) {
        final List<StructuralVariantAnnotation> annotations = annotator.annotateVariants(variants);

        final List<StructuralVariantAnnotation> copy = Lists.newArrayList(annotations);
        final List<GeneFusion> fusions = processFusions(copy);
        final List<GeneDisruption> disruptions = processDisruptions(copy);

        return ImmutableStructuralVariantAnalysis.of(annotations, fusions, disruptions);
    }

    @NotNull
    private List<GeneFusion> processFusions(final List<StructuralVariantAnnotation> annotations) {
        // NERA: left is upstream, right is downstream
        final List<List<Pair<Transcript, Transcript>>> fusionsPerVariant = Lists.newArrayList();
        for (final StructuralVariantAnnotation annotation : annotations) {
            final List<Pair<Transcript, Transcript>> fusions = Lists.newArrayList();

            for (final GeneAnnotation startGene : annotation.start()) {
                final boolean startUpstream = isUpstream(startGene);

                for (final GeneAnnotation endGene : annotation.end()) {
                    final boolean endUpstream = isUpstream(endGene);
                    if (startUpstream == endUpstream) {
                        continue;
                    }

                    for (final Transcript t1 : intronic(startGene.transcripts())) {
                        for (final Transcript t2 : intronic(endGene.transcripts())) {
                            if (!isPotentiallyRelevantFusion(t1, t2)) {
                                continue;
                            }

                            if (startUpstream && t1.exonUpstreamPhase() == t2.exonDownstreamPhase()) {
                                fusions.add(Pair.of(t1, t2));
                            } else if (!startUpstream && t2.exonUpstreamPhase() == t1.exonDownstreamPhase()) {
                                fusions.add(Pair.of(t2, t1));
                            }
                        }
                    }
                }
            }

            fusionsPerVariant.add(fusions);
        }

        return toReportableGeneFusions(fusionsPerVariant);
    }

    private static boolean isPotentiallyRelevantFusion(@NotNull Transcript t1, @NotNull Transcript t2) {
        final boolean sameGene = t1.geneName().equals(t2.geneName());
        if (sameGene) {
            // NERA: skip fusions between different transcripts in the same gene,
            if (!t1.transcriptId().equals(t2.transcriptId())) {
                return false;
            }
            // NERA: skip fusions within the same intron
            if (intronicDisruption(t1, t2)) {
                return false;
            }
        }
        return true;
    }

    @NotNull
    private List<GeneFusion> toReportableGeneFusions(@NotNull List<List<Pair<Transcript, Transcript>>> fusionsPerVariant) {
        final List<GeneFusion> result = Lists.newArrayList();
        for (final List<Pair<Transcript, Transcript>> fusions : fusionsPerVariant) {
            Optional<Pair<Transcript, Transcript>> reportableFusion = determineReportableFusion(fusions);

            for (final Pair<Transcript, Transcript> fusion : fusions) {
                final Transcript upstream = fusion.getLeft();
                final Transcript downstream = fusion.getRight();

                final CosmicFusionData cosmic = transcriptsMatchKnownFusion(upstream, downstream);
                final boolean promiscuousEnd = oneEndPromiscuous(upstream, downstream);
                final boolean reportable =
                        reportableFusion.isPresent() && reportableFusion.get() == fusion && (cosmic != null || promiscuousEnd);

                final GeneFusion geneFusion = ImmutableGeneFusion.builder()
                        .reportable(reportable)
                        .upstreamLinkedAnnotation(upstream)
                        .downstreamLinkedAnnotation(downstream)
                        .cosmicURL(cosmic != null ? cosmic.cosmicURL() : "")
                        .build();

                result.add(geneFusion);
            }
        }
        return result;
    }

    @NotNull
    private static Optional<Pair<Transcript, Transcript>> determineReportableFusion(@NotNull List<Pair<Transcript, Transcript>> fusions) {
        // NERA: Select either the canonical -> canonical transcript fusion
        //  then the one with the most exons where one end is canonical
        //  then the one with the most exons combined transcript

        Optional<Pair<Transcript, Transcript>> reportableFusion =
                fusions.stream().filter(pair -> pair.getLeft().isCanonical() && pair.getRight().isCanonical()).findFirst();

        if (!reportableFusion.isPresent()) {
            reportableFusion = fusions.stream()
                    .filter(pair -> pair.getLeft().isCanonical() || pair.getRight().isCanonical())
                    .sorted(Comparator.comparingInt(a -> a.getLeft().exonMax() + a.getRight().exonMax()))
                    .reduce((a, b) -> b);
        }

        if (!reportableFusion.isPresent()) {
            reportableFusion = fusions.stream()
                    .sorted(Comparator.comparingInt(a -> a.getLeft().exonMax() + a.getRight().exonMax()))
                    .reduce((a, b) -> b);
        }

        return reportableFusion;
    }

    @NotNull
    private List<GeneDisruption> processDisruptions(final List<StructuralVariantAnnotation> annotations) {
        final List<GeneAnnotation> geneAnnotations = Lists.newArrayList();
        for (final StructuralVariantAnnotation annotation : annotations) {
            @SuppressWarnings("ConstantConditions")
            final boolean intronicExists = annotation.start()
                    .stream()
                    .filter(gene -> gene.canonical() != null)
                    .anyMatch(gene -> annotation.end()
                            .stream()
                            .filter(o -> o.canonical() != null)
                            .anyMatch(o -> intronicDisruption(gene.canonical(), o.canonical())));
            if (intronicExists && annotation.variant().type() != StructuralVariantType.INV) {
                continue;
            }

            geneAnnotations.addAll(annotation.annotations());
        }

        final Multimap<String, GeneAnnotation> geneMap = ArrayListMultimap.create();
        geneAnnotations.forEach(g -> geneMap.put(g.geneName(), g));

        final List<GeneDisruption> disruptions = Lists.newArrayList();
        for (final String geneName : geneMap.keySet()) {
            for (final GeneAnnotation g : geneMap.get(geneName)) {
                for (final Transcript transcript : g.transcripts()) {
                    final GeneDisruption disruption = ImmutableGeneDisruption.builder()
                            .reportable(inHmfPanel(g) && transcript.isCanonical())
                            .linkedAnnotation(transcript)
                            .build();

                    disruptions.add(disruption);
                }
            }
        }

        return disruptions;
    }

    private boolean inHmfPanel(@NotNull GeneAnnotation gene) {
        return hmfGenePanelRegions.stream().anyMatch(region -> gene.synonyms().contains(region.geneID()));
    }

    private boolean transcriptsMatchKnownFusion(@NotNull CosmicFusionData fusion, @NotNull Transcript five, @NotNull Transcript three) {
        String fiveTranscript = fusion.fiveTranscript();
        String threeTranscript = fusion.threeTranscript();

        final boolean fiveValid = fiveTranscript == null
                ? five.parent().synonyms().stream().anyMatch(s -> s.equals(fusion.fiveGene()))
                : fiveTranscript.equals(five.transcriptId());
        final boolean threeValid = threeTranscript == null
                ? three.parent().synonyms().stream().anyMatch(s -> s.equals(fusion.threeGene()))
                : threeTranscript.equals(three.transcriptId());
        return fiveValid && threeValid;
    }

    @Nullable
    private CosmicFusionData transcriptsMatchKnownFusion(@NotNull Transcript five, @NotNull Transcript three) {
        return cosmicFusionModel.fusions().stream().filter(f -> transcriptsMatchKnownFusion(f, five, three)).findFirst().orElse(null);
    }

    private boolean oneEndPromiscuous(@NotNull Transcript five, @NotNull Transcript three) {
        @SuppressWarnings("ConstantConditions")
        final boolean promiscuousFive = cosmicFusionModel.promiscuousFivePrime()
                .stream()
                .anyMatch(p -> p.transcript() != null ? p.transcript().equals(five.transcriptId()) : p.geneName().equals(five.geneName()));
        @SuppressWarnings("ConstantConditions")
        final boolean promiscuousThree = cosmicFusionModel.promiscuousThreePrime()
                .stream()
                .anyMatch(p -> p.transcript() != null
                        ? p.transcript().equals(three.transcriptId())
                        : p.geneName().equals(three.geneName()));
        return promiscuousFive || promiscuousThree;
    }

    private static boolean intronicDisruption(@NotNull Transcript t1, @NotNull Transcript t2) {
        final boolean sameTranscript = t1.transcriptId().equals(t2.transcriptId());
        final boolean bothIntronic = t1.isIntronic() && t2.isIntronic();
        final boolean sameExonUpstream = t1.exonUpstream() == t2.exonUpstream();
        return sameTranscript && bothIntronic && sameExonUpstream;
    }

    private static boolean isUpstream(@NotNull GeneAnnotation annotation) {
        int orientation = annotation.variant().orientation(annotation.isStart());
        return annotation.strand() * orientation > 0;
    }

    @NotNull
    private static List<Transcript> intronic(@NotNull List<Transcript> transcripts) {
        return transcripts.stream().filter(Transcript::isIntronic).collect(Collectors.toList());
    }
}
