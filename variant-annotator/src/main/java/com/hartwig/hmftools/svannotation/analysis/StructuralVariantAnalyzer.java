package com.hartwig.hmftools.svannotation.analysis;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class StructuralVariantAnalyzer {
    private static final int EXON_THRESHOLD = 1;

    @NotNull
    private final VariantAnnotator annotator;
    @NotNull
    private final Collection<HmfGenomeRegion> hmfGenePanelRegions;
    @NotNull
    private final KnownFusionsModel knownFusionsModel;

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnalyzer.class);

    public StructuralVariantAnalyzer(@NotNull final VariantAnnotator annotator,
            @NotNull final Collection<HmfGenomeRegion> hmfGenePanelRegions, @NotNull final KnownFusionsModel knownFusionsModel) {
        this.annotator = annotator;
        this.hmfGenePanelRegions = hmfGenePanelRegions;
        this.knownFusionsModel = knownFusionsModel;
    }

    @NotNull
    public StructuralVariantAnalysis run(@NotNull final List<EnrichedStructuralVariant> variants) {
        LOGGER.debug("annotating variants");
        final List<StructuralVariantAnnotation> annotations = annotator.annotateVariants(variants);

        final List<StructuralVariantAnnotation> copy = Lists.newArrayList(annotations);

        LOGGER.debug("calling process fusions");
        final List<GeneFusion> fusions = processFusions(copy);

        LOGGER.debug("calling process disruptions");
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
            if (intronicDisruptionOnSameTranscript(t1, t2)) {
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
                final boolean matchesKnownFusion = transcriptsMatchKnownFusion(knownFusionsModel, upstream, downstream);
                final boolean reportable = reportableFusion.isPresent() && reportableFusion.get() == fusion && matchesKnownFusion && (
                        (postCoding(downstream) != null && !postCoding(downstream)) && (postCoding(upstream) == null
                                || !postCoding(upstream)) && !(intragenic(upstream, downstream) && upstream.exonUpstreamPhase() == -1));
                final GeneFusion geneFusion = ImmutableGeneFusion.builder()
                        .reportable(reportable)
                        .upstreamLinkedAnnotation(upstream)
                        .downstreamLinkedAnnotation(downstream)
                        .primarySource(knownFusionsModel.primarySource(upstream.parent().synonyms(), downstream.parent().synonyms()))
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
            final boolean pureIntronicDisruptionCanonical = annotation.start()
                    .stream()
                    .filter(gene -> gene.canonical() != null)
                    .anyMatch(gene -> annotation.end()
                            .stream()
                            .filter(other -> other.canonical() != null)
                            .anyMatch(other -> intronicDisruptionOnSameTranscript(gene.canonical(), other.canonical())));
            if (pureIntronicDisruptionCanonical && annotation.variant().type() != StructuralVariantType.INV) {
                continue;
            }

            geneAnnotations.addAll(annotation.annotations());
        }

        final Multimap<String, GeneAnnotation> geneMap = ArrayListMultimap.create();
        geneAnnotations.forEach(g -> geneMap.put(g.geneName(), g));

        final List<GeneDisruption> disruptions = Lists.newArrayList();
        for (final String geneName : geneMap.keySet()) {
            for (final GeneAnnotation gene : geneMap.get(geneName)) {
                for (final Transcript transcript : gene.transcripts()) {
                    final GeneDisruption disruption = ImmutableGeneDisruption.builder()
                            .reportable(inHmfPanel(gene) && transcript.isCanonical())
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

    private static boolean transcriptsMatchKnownFusion(@NotNull final KnownFusionsModel fusionsModel, @NotNull final Transcript five,
            @NotNull final Transcript three) {
        return fusionsModel.exactMatch(five.parent().synonyms(), three.parent().synonyms())
                || fusionsModel.intergenicPromiscuousMatch(five.parent().synonyms(), three.parent().synonyms()) || (
                fusionsModel.intragenicPromiscuousMatch(five.parent().synonyms(), three.parent().synonyms())
                        && three.exonDownstream() - five.exonUpstream() > EXON_THRESHOLD);
    }

    private static boolean intronicDisruptionOnSameTranscript(@NotNull Transcript t1, @NotNull Transcript t2) {
        final boolean sameTranscript = t1.transcriptId().equals(t2.transcriptId());
        final boolean bothIntronic = t1.isIntronic() && t2.isIntronic();
        final boolean sameExonUpstream = t1.exonUpstream() == t2.exonUpstream();
        return sameTranscript && bothIntronic && sameExonUpstream;
    }

    private static boolean isUpstream(@NotNull GeneAnnotation gene) {
        int orientation = gene.variant().orientation(gene.isStart());
        return gene.strand() * orientation > 0;
    }

    @NotNull
    private static List<Transcript> intronic(@NotNull List<Transcript> transcripts) {
        return transcripts.stream().filter(Transcript::isIntronic).collect(Collectors.toList());
    }

    private static boolean preCoding(@NotNull final Transcript transcript) {
        final int strand = transcript.parent().strand();
        final long position = transcript.parent().variant().position(transcript.parent().isStart());
        return (strand == 1 && position < transcript.codingStart()) || (strand == -1 && position > transcript.codingEnd());
    }

    @Nullable
    private static Boolean postCoding(@NotNull final Transcript transcript) {
        if (transcript.codingEnd() == null || transcript.codingStart() == null) {
            return null;
        }
        final int strand = transcript.parent().strand();
        final long position = transcript.parent().variant().position(transcript.parent().isStart());
        return (strand == 1 && position > transcript.codingEnd()) || (strand == -1 && position < transcript.codingStart());
    }

    private static boolean intragenic(@NotNull final Transcript upstream, @NotNull final Transcript downstream) {
        return upstream.parent().synonyms().stream().anyMatch(downstream.parent().synonyms()::contains);
    }
}
