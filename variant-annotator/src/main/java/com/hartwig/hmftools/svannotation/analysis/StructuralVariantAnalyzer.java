package com.hartwig.hmftools.svannotation.analysis;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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

public class StructuralVariantAnalyzer {

    private final VariantAnnotator annotator;
    private final Collection<HmfGenomeRegion> regions;
    private final COSMICGeneFusionModel fusionModel;

    public StructuralVariantAnalyzer(final VariantAnnotator annotator, final Collection<HmfGenomeRegion> regions,
            final COSMICGeneFusionModel fusionModel) {
        this.annotator = annotator;
        this.regions = regions;
        this.fusionModel = fusionModel;
    }

    private boolean inHmfPanel(final GeneAnnotation g) {
        return regions.stream().anyMatch(r -> g.getSynonyms().contains(r.geneID()));
    }

    private boolean transcriptsMatchKnownFusion(final COSMICGeneFusionData fusion, final Transcript five, final Transcript three) {
        final boolean fiveValid = fusion.fiveTranscript() == null
                ? five.getGeneAnnotation().getSynonyms().stream().anyMatch(s -> s.equals(fusion.fiveGene()))
                : fusion.fiveTranscript().equals(five.getTranscriptId());
        final boolean threeValid = fusion.threeTranscript() == null ? three.getGeneAnnotation()
                .getSynonyms()
                .stream()
                .anyMatch(s -> s.equals(fusion.threeGene())) : fusion.threeTranscript().equals(three.getTranscriptId());
        return fiveValid && threeValid;
    }

    private COSMICGeneFusionData transcriptsMatchKnownFusion(final Transcript five, final Transcript three) {
        return fusionModel.fusions().stream().filter(f -> transcriptsMatchKnownFusion(f, five, three)).findFirst().orElse(null);
    }

    private boolean isPromiscuous(final GeneAnnotation gene) {
        return Stream.of(fusionModel.promiscuousFivePrime(), fusionModel.promiscuousThreePrime())
                .anyMatch(l -> l.stream().anyMatch(g -> gene.getSynonyms().contains(g.geneName())));
    }

    private boolean oneEndPromiscuous(final Transcript five, final Transcript three) {
        final boolean promiscuousFive = fusionModel.promiscuousFivePrime()
                .stream()
                .anyMatch(p -> p.transcript() != null
                        ? p.transcript().equals(five.getTranscriptId())
                        : p.geneName().equals(five.getGeneName()));
        final boolean promiscuousThree = fusionModel.promiscuousThreePrime()
                .stream()
                .anyMatch(p -> p.transcript() != null
                        ? p.transcript().equals(three.getTranscriptId())
                        : p.geneName().equals(three.getGeneName()));
        return promiscuousFive || promiscuousThree;
    }

    private boolean intronicDisruption(final Transcript a, final Transcript b) {
        final boolean sameTranscript = a.getTranscriptId().equals(b.getTranscriptId());
        final boolean bothIntronic = a.isIntronic() && b.isIntronic();
        final boolean sameExonUpstream = a.getExonUpstream() == b.getExonUpstream();
        return sameTranscript && bothIntronic && sameExonUpstream;
    }

    private static int orientation(final GeneAnnotation g) {
        final StructuralVariant sv = g.getVariant();
        return sv.orientation(g.isStart());
    }

    private static List<Transcript> intronic(final List<Transcript> list) {
        return list.stream().filter(Transcript::isIntronic).collect(Collectors.toList());
    }

    private List<GeneFusion> processFusions(final List<StructuralVariantAnnotation> annotations) {

        // left is upstream, right is downstream

        final List<List<Pair<Transcript, Transcript>>> fusionsPerVariant = Lists.newArrayList();
        for (final StructuralVariantAnnotation sv : annotations) {

            final List<Pair<Transcript, Transcript>> fusions = Lists.newArrayList();

            for (final GeneAnnotation g : sv.getStart()) {

                final boolean g_upstream = g.getStrand() * orientation(g) > 0;

                for (final GeneAnnotation o : sv.getEnd()) {

                    final boolean o_upstream = o.getStrand() * orientation(o) > 0;
                    if (g_upstream == o_upstream) {
                        continue;
                    }

                    for (final Transcript t1 : intronic(g.getTranscripts())) {
                        for (final Transcript t2 : intronic(o.getTranscripts())) {

                            final boolean sameGene = t1.getGeneName().equals(t2.getGeneName());
                            if (sameGene) {
                                // skip fusions between different transcripts in the same gene,
                                if (!t1.getTranscriptId().equals(t2.getTranscriptId())) {
                                    continue;
                                }
                                // skip fusions within the same intron
                                if (intronicDisruption(t1, t2)) {
                                    continue;
                                }
                            }

                            if (g_upstream && t1.getExonUpstreamPhase() == t2.getExonDownstreamPhase()) {
                                fusions.add(Pair.of(t1, t2));
                            } else if (!g_upstream && t2.getExonUpstreamPhase() == t1.getExonDownstreamPhase()) {
                                fusions.add(Pair.of(t2, t1));
                            }

                        }
                    }
                }
            }

            fusionsPerVariant.add(fusions);
        }

        // transform results to reported details

        final List<GeneFusion> result = Lists.newArrayList();
        for (final List<Pair<Transcript, Transcript>> fusions : fusionsPerVariant) {

            // from here, select either the canonical -> canonical transcript fusion
            // then the longest where one end is canonical
            // then the longest combined transcript

            Optional<Pair<Transcript, Transcript>> reportableFusion =
                    fusions.stream().filter(p -> p.getLeft().isCanonical() && p.getRight().isCanonical()).findFirst();

            if (!reportableFusion.isPresent()) {
                reportableFusion = fusions.stream()
                        .filter(p -> p.getLeft().isCanonical() || p.getRight().isCanonical())
                        .sorted(Comparator.comparingInt(a -> a.getLeft().getExonMax() + a.getRight().getExonMax()))
                        .reduce((a, b) -> b); // get longest
            }

            if (!reportableFusion.isPresent()) {
                reportableFusion = fusions.stream()
                        .sorted(Comparator.comparingInt(a -> a.getLeft().getExonMax() + a.getRight().getExonMax()))
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

    private List<GeneDisruption> processDisruptions(final List<StructuralVariantAnnotation> annotations) {

        final List<GeneAnnotation> geneAnnotations = Lists.newArrayList();
        for (final StructuralVariantAnnotation sv : annotations) {

            final boolean intronicExists = sv.getStart()
                    .stream()
                    .filter(g -> g.getCanonical() != null)
                    .anyMatch(g -> sv.getEnd()
                            .stream()
                            .filter(o -> o.getCanonical() != null)
                            .anyMatch(o -> intronicDisruption(g.getCanonical(), o.getCanonical())));
            if (intronicExists && sv.getVariant().type() != StructuralVariantType.INV) {
                continue;
            }

            geneAnnotations.addAll(sv.getAnnotations());
        }

        final ArrayListMultimap<String, GeneAnnotation> geneMap = ArrayListMultimap.create();
        geneAnnotations.forEach(g -> geneMap.put(g.getGeneName(), g));

        final List<GeneDisruption> disruptions = Lists.newArrayList();
        for (final String geneName : geneMap.keySet()) {
            for (final GeneAnnotation g : geneMap.get(geneName)) {

                for (final Transcript transcript : g.getTranscripts()) {

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
