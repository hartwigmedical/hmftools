package com.hartwig.hmftools.svannotation.analysis;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cosmicfusions.COSMICGeneFusionData;
import com.hartwig.hmftools.common.cosmicfusions.COSMICGeneFusionModel;
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
    private final Collection<HmfGenomeRegion> regions;
    @NotNull
    private final COSMICGeneFusionModel fusionModel;

    public StructuralVariantAnalyzer(@NotNull final VariantAnnotator annotator, @NotNull final Collection<HmfGenomeRegion> regions,
            @NotNull final COSMICGeneFusionModel fusionModel) {
        this.annotator = annotator;
        this.regions = regions;
        this.fusionModel = fusionModel;
    }

    private boolean inHmfPanel(final GeneAnnotation gene) {
        return regions.stream().anyMatch(region -> gene.synonyms().contains(region.geneID()));
    }

    private boolean transcriptsMatchKnownFusion(final COSMICGeneFusionData fusion, final Transcript five, final Transcript three) {
        String fiveTranscript = fusion.fiveTranscript();
        String threeTranscript = fusion.threeTranscript();

        final boolean fiveValid = fiveTranscript == null
                ? five.geneAnnotation().synonyms().stream().anyMatch(s -> s.equals(fusion.fiveGene()))
                : fiveTranscript.equals(five.transcriptId());
        final boolean threeValid = threeTranscript == null
                ? three.geneAnnotation().synonyms().stream().anyMatch(s -> s.equals(fusion.threeGene()))
                : threeTranscript.equals(three.transcriptId());
        return fiveValid && threeValid;
    }

    @Nullable
    private COSMICGeneFusionData transcriptsMatchKnownFusion(final Transcript five, final Transcript three) {
        return fusionModel.fusions().stream().filter(f -> transcriptsMatchKnownFusion(f, five, three)).findFirst().orElse(null);
    }

    private boolean oneEndPromiscuous(final Transcript five, final Transcript three) {
        @SuppressWarnings("ConstantConditions")
        final boolean promiscuousFive = fusionModel.promiscuousFivePrime()
                .stream()
                .anyMatch(p -> p.transcript() != null ? p.transcript().equals(five.transcriptId()) : p.geneName().equals(five.geneName()));
        @SuppressWarnings("ConstantConditions")
        final boolean promiscuousThree = fusionModel.promiscuousThreePrime()
                .stream()
                .anyMatch(p -> p.transcript() != null
                        ? p.transcript().equals(three.transcriptId())
                        : p.geneName().equals(three.geneName()));
        return promiscuousFive || promiscuousThree;
    }

    private boolean intronicDisruption(final Transcript a, final Transcript b) {
        final boolean sameTranscript = a.transcriptId().equals(b.transcriptId());
        final boolean bothIntronic = a.isIntronic() && b.isIntronic();
        final boolean sameExonUpstream = a.exonUpstream() == b.exonUpstream();
        return sameTranscript && bothIntronic && sameExonUpstream;
    }

    private static int orientation(final GeneAnnotation g) {
        final StructuralVariant sv = g.variant();
        return sv.orientation(g.isStart());
    }

    @NotNull
    private static List<Transcript> intronic(final List<Transcript> list) {
        return list.stream().filter(Transcript::isIntronic).collect(Collectors.toList());
    }

    private List<GeneFusion> processFusions(final List<StructuralVariantAnnotation> annotations) {
        // NERA: left is upstream, right is downstream
        final List<List<Pair<Transcript, Transcript>>> fusionsPerVariant = Lists.newArrayList();
        for (final StructuralVariantAnnotation sv : annotations) {

            final List<Pair<Transcript, Transcript>> fusions = Lists.newArrayList();

            for (final GeneAnnotation g : sv.start()) {

                final boolean g_upstream = g.strand() * orientation(g) > 0;

                for (final GeneAnnotation o : sv.end()) {

                    final boolean o_upstream = o.strand() * orientation(o) > 0;
                    if (g_upstream == o_upstream) {
                        continue;
                    }

                    for (final Transcript t1 : intronic(g.transcripts())) {
                        for (final Transcript t2 : intronic(o.transcripts())) {

                            final boolean sameGene = t1.geneName().equals(t2.geneName());
                            if (sameGene) {
                                // NERA: skip fusions between different transcripts in the same gene,
                                if (!t1.transcriptId().equals(t2.transcriptId())) {
                                    continue;
                                }
                                // NERA: skip fusions within the same intron
                                if (intronicDisruption(t1, t2)) {
                                    continue;
                                }
                            }

                            if (g_upstream && t1.exonUpstreamPhase() == t2.exonDownstreamPhase()) {
                                fusions.add(Pair.of(t1, t2));
                            } else if (!g_upstream && t2.exonUpstreamPhase() == t1.exonDownstreamPhase()) {
                                fusions.add(Pair.of(t2, t1));
                            }

                        }
                    }
                }
            }

            fusionsPerVariant.add(fusions);
        }

        // NERA: transform results to reported details

        final List<GeneFusion> result = Lists.newArrayList();
        for (final List<Pair<Transcript, Transcript>> fusions : fusionsPerVariant) {

            // NERA: from here, select either the canonical -> canonical transcript fusion
            // then the longest where one end is canonical
            // then the longest combined transcript

            Optional<Pair<Transcript, Transcript>> reportableFusion =
                    fusions.stream().filter(p -> p.getLeft().isCanonical() && p.getRight().isCanonical()).findFirst();

            if (!reportableFusion.isPresent()) {
                reportableFusion = fusions.stream()
                        .filter(p -> p.getLeft().isCanonical() || p.getRight().isCanonical())
                        .sorted(Comparator.comparingInt(a -> a.getLeft().exonMax() + a.getRight().exonMax()))
                        .reduce((a, b) -> b); // get longest
            }

            if (!reportableFusion.isPresent()) {
                reportableFusion = fusions.stream().sorted(Comparator.comparingInt(a -> a.getLeft().exonMax() + a.getRight().exonMax()))
                        .reduce((a, b) -> b); // get longest
            }

            for (final Pair<Transcript, Transcript> fusion : fusions) {
                final Transcript upstream = fusion.getLeft(), downstream = fusion.getRight();

                final COSMICGeneFusionData cosmic = transcriptsMatchKnownFusion(upstream, downstream);
                final boolean promiscuousEnd = oneEndPromiscuous(upstream, downstream);
                final boolean reportable =
                        reportableFusion.isPresent() && reportableFusion.get() == fusion && (cosmic != null || promiscuousEnd);

                final GeneFusion details = ImmutableGeneFusion.builder()
                        .reportable(reportable)
                        .upstreamLinkedAnnotation(upstream)
                        .downstreamLinkedAnnotation(downstream)
                        .cosmicURL(cosmic != null ? cosmic.cosmicURL() : "")
                        .build();

                result.add(details);
            }

        }

        return result;
    }

    @NotNull
    private List<GeneDisruption> processDisruptions(final List<StructuralVariantAnnotation> annotations) {
        final List<GeneAnnotation> geneAnnotations = Lists.newArrayList();
        for (final StructuralVariantAnnotation sv : annotations) {

            final boolean intronicExists = sv.start()
                    .stream()
                    .filter(g -> g.canonical() != null)
                    .anyMatch(g -> sv.end()
                            .stream()
                            .filter(o -> o.canonical() != null)
                            .anyMatch(o -> intronicDisruption(g.canonical(), o.canonical())));
            if (intronicExists && sv.variant().type() != StructuralVariantType.INV) {
                continue;
            }

            geneAnnotations.addAll(sv.annotations());
        }

        final ArrayListMultimap<String, GeneAnnotation> geneMap = ArrayListMultimap.create();
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

    @NotNull
    public StructuralVariantAnalysis run(@NotNull final List<StructuralVariant> variants) {
        final List<StructuralVariantAnnotation> annotations = annotator.annotateVariants(variants);

        final List<StructuralVariantAnnotation> copy = Lists.newArrayList(annotations);
        final List<GeneFusion> fusions = processFusions(copy);
        final List<GeneDisruption> disruptions = processDisruptions(copy);

        return ImmutableStructuralVariantAnalysis.of(annotations, fusions, disruptions);
    }
}
