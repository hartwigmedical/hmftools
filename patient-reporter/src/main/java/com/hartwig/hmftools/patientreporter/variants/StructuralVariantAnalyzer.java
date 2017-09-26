package com.hartwig.hmftools.patientreporter.variants;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.patientreporter.data.COSMICGeneFusionData;
import com.hartwig.hmftools.patientreporter.data.COSMICGeneFusionModel;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.svannotation.GeneAnnotation;
import com.hartwig.hmftools.svannotation.Transcript;
import com.hartwig.hmftools.svannotation.VariantAnnotation;
import com.hartwig.hmftools.svannotation.VariantAnnotator;

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

    private boolean transcriptsMatchKnownFusion(final Transcript five, final Transcript three) {
        return fusionModel.fusions().stream().anyMatch(f -> transcriptsMatchKnownFusion(f, five, three));
    }

    private boolean isPromiscuous(final GeneAnnotation gene) {
        return Stream.of(fusionModel.promiscuousFivePrime(), fusionModel.promiscuousThreePrime())
                .anyMatch(l -> l.stream().anyMatch(g -> gene.getSynonyms().contains(g.GeneName())));
    }

    private boolean oneEndPromiscuous(final Transcript five, final Transcript three) {
        final boolean promiscuousFive = fusionModel.promiscuousFivePrime()
                .stream()
                .anyMatch(p -> p.Transcript() != null
                        ? p.Transcript().equals(five.getTranscriptId())
                        : p.GeneName().equals(five.getGeneName()));
        final boolean promiscuousThree = fusionModel.promiscuousThreePrime()
                .stream()
                .anyMatch(p -> p.Transcript() != null
                        ? p.Transcript().equals(three.getTranscriptId())
                        : p.GeneName().equals(three.getGeneName()));
        return promiscuousFive || promiscuousThree;
    }

    private boolean intronicDisruption(final Transcript a, final Transcript b) {
        final boolean sameTranscript = a.getTranscriptId().equals(b.getTranscriptId());
        final boolean bothIntronic = a.isIntronic() && b.isIntronic();
        final boolean sameExonUpstream = a.getExonUpstream() == b.getExonUpstream();
        return sameTranscript && bothIntronic && sameExonUpstream;
    }

    private String exonDescription(final Transcript t) {
        if (t.isPromoter()) {
            return "Promoter Region";
        } else if (t.isExonic()) {
            return String.format("Exon %d", t.getExonUpstream());
        } else if (t.isIntronic()) {
            return String.format("Intron %d", t.getExonUpstream());
        } else {
            return String.format("Error up(%d) down(%d)", t.getExonUpstream(), t.getExonDownstream());
        }
    }

    private String exonSelection(final Transcript t, final boolean upstream) {
        return String.format("Exon %d", upstream ? t.getExonUpstream() : t.getExonDownstream());
    }

    private List<StructuralVariantAnalysis.GeneFusion> processFusions(final List<VariantAnnotation> annotations) {

        // left is upstream, right is downstream
        final List<Pair<Transcript, Transcript>> fusions = Lists.newArrayList();


        for (final VariantAnnotation sv : annotations) {

            final List<Pair<Transcript, Transcript>> svFusions = Lists.newArrayList();

            for (final GeneAnnotation g : sv.getStart().getGeneAnnotations()) {

                final boolean g_upstream = g.getStrand() * g.getBreakend().getOrientation() > 0;

                for (final GeneAnnotation o : sv.getEnd().getGeneAnnotations()) {

                    final boolean o_upstream = o.getStrand() * o.getBreakend().getOrientation() > 0;
                    if (g_upstream == o_upstream) {
                        continue;
                    }

                    for (final Transcript t1 : g.getTranscripts()) {
                        if (!t1.isIntronic()) {
                            continue;
                        }

                        for (final Transcript t2 : o.getTranscripts()) {
                            if (!t2.isIntronic()) {
                                continue;
                            }

                            if (g_upstream && t1.getExonUpstreamPhase() == t2.getExonDownstreamPhase()) {
                                svFusions.add(Pair.of(t1, t2));
                            } else if (!g_upstream && t2.getExonUpstreamPhase() == t1.getExonDownstreamPhase()) {
                                svFusions.add(Pair.of(t2, t1));
                            }

                        }
                    }
                }
            }

            // from here, select either the canonical -> canonical transcript fusion
            // then the longest where one end is canonical
            // then the longest combined transcript

            Optional<Pair<Transcript, Transcript>> fusion =
                    svFusions.stream().filter(p -> p.getLeft().isCanonical() && p.getRight().isCanonical()).findFirst();

            if (fusion.isPresent()) {
                fusions.add(fusion.get());
                continue;
            }

            fusion = svFusions.stream()
                    .filter(p -> p.getLeft().isCanonical() || p.getRight().isCanonical())
                    .sorted(Comparator.comparingInt(a -> a.getLeft().getExonMax() + a.getRight().getExonMax()))
                    .reduce((a, b) -> b); // get longest

            if (fusion.isPresent()) {
                fusions.add(fusion.get());
                continue;
            }

            svFusions.stream()
                    .sorted(Comparator.comparingInt(a -> a.getLeft().getExonMax() + a.getRight().getExonMax()))
                    .reduce((a, b) -> b) // get longest
                    .ifPresent(fusions::add);
        }

        // transform results to reported details

        final List<StructuralVariantAnalysis.GeneFusion> result = Lists.newArrayList();
        for (final Pair<Transcript, Transcript> fusion : fusions) {
            final Transcript upstream = fusion.getLeft(), downstream = fusion.getRight();
            final boolean sameGene = upstream.getGeneName().equals(downstream.getGeneName());

            if (sameGene) {
                if (!intronicDisruption(upstream, downstream) && isPromiscuous(upstream.getGeneAnnotation())) {
                    // okay
                } else {
                    continue;
                }
            } else if (transcriptsMatchKnownFusion(upstream, downstream)) {
                // in cosmic fusion list
            } else if (oneEndPromiscuous(upstream, downstream)) {
                // one end is promiscuous
            } else {
                continue;
            }

            final Double fiveAF = upstream.getBreakend().getAlleleFrequency();
            final Double threeAF = downstream.getBreakend().getAlleleFrequency();

            final StructuralVariantAnalysis.GeneFusion details = ImmutableGeneFusion.builder()
                    .Type(upstream.getBreakend().getStructuralVariant().getVariant().type().toString())
                    .Start(upstream.getBreakend().getPositionString())
                    .GeneStart(upstream.getGeneName())
                    .GeneContextStart(exonSelection(upstream, true))
                    .TranscriptStart(upstream.getTranscriptId())
                    .End(downstream.getBreakend().getPositionString())
                    .GeneEnd(downstream.getGeneName())
                    .GeneContextEnd(exonSelection(downstream, false))
                    .TranscriptEnd(downstream.getTranscriptId())
                    .VAF(PatientReportFormat.formatNullablePercent(fiveAF) + " " + PatientReportFormat.formatNullablePercent(threeAF))
                    .build();

            result.add(details);
            annotations.remove(upstream.getBreakend().getStructuralVariant());
        }

        return result;
    }

    private List<StructuralVariantAnalysis.GeneDisruption> processDisruptions(final List<VariantAnnotation> annotations) {

        final List<GeneAnnotation> geneAnnotations = Lists.newArrayList();
        for (final VariantAnnotation sv : annotations) {

            final boolean intronicExists = sv.getStart()
                    .getGeneAnnotations()
                    .stream()
                    .filter(g -> g.getCanonical() != null)
                    .anyMatch(g -> sv.getEnd()
                            .getGeneAnnotations()
                            .stream()
                            .filter(o -> o.getCanonical() != null)
                            .anyMatch(o -> intronicDisruption(g.getCanonical(), o.getCanonical())));
            if (intronicExists) {
                continue;
            }

            geneAnnotations.addAll(sv.getStart().getGeneAnnotations());
            geneAnnotations.addAll(sv.getEnd().getGeneAnnotations());
        }

        final ArrayListMultimap<String, GeneAnnotation> geneMap = ArrayListMultimap.create();
        for (final GeneAnnotation g : geneAnnotations) {
            if (!inHmfPanel(g)) {
                continue;
            }
            geneMap.put(g.getGeneName(), g);
        }

        final List<StructuralVariantAnalysis.GeneDisruption> disruptions = Lists.newArrayList();
        for (final String geneName : geneMap.keySet()) {
            for (final GeneAnnotation g : geneMap.get(geneName)) {

                // don't care if we aren't in the canonical transcript
                if (g.getCanonical() == null) {
                    continue;
                }

                final StructuralVariantAnalysis.GeneDisruption disruption = ImmutableGeneDisruption.builder()
                        .GeneName(geneName)
                        .Location(g.getBreakend().getPositionString())
                        .GeneContext(exonDescription(g.getCanonical()))
                        .Transcript(g.getCanonical().getTranscriptId())
                        .Partner(g.getOtherBreakend().getPositionString())
                        .Type(g.getBreakend().getStructuralVariant().getVariant().type().toString())
                        .Orientation(g.getBreakend().getOrientation() > 0 ? "5'" : "3'")
                        .VAF(PatientReportFormat.formatNullablePercent(g.getBreakend().getAlleleFrequency()))
                        .build();

                disruptions.add(disruption);
                annotations.remove(g.getBreakend().getStructuralVariant());
            }
        }

        return disruptions;
    }

    @NotNull
    public StructuralVariantAnalysis run(@NotNull final List<StructuralVariant> variants) {

        final List<VariantAnnotation> annotations = annotator.annotateVariants(variants);
        final List<StructuralVariantAnalysis.GeneFusion> fusions = processFusions(annotations);
        final List<StructuralVariantAnalysis.GeneDisruption> disruptions = processDisruptions(annotations);

        return ImmutableStructuralVariantAnalysis.of(annotations, fusions, disruptions);
    }
}
