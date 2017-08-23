package com.hartwig.hmftools.patientreporter.variants;

import java.util.Comparator;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.slicing.HmfSlicer;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.svannotation.GeneAnnotation;
import com.hartwig.hmftools.svannotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.svannotation.StructuralVariantAnnotator;
import com.hartwig.hmftools.svannotation.TranscriptAnnotation;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

public class StructuralVariantAnalyzer {

    private final StructuralVariantAnnotator annotator;
    private final HmfSlicer slicer;

    public StructuralVariantAnalyzer(final StructuralVariantAnnotator annotator, final HmfSlicer slicer) {
        this.annotator = annotator;
        this.slicer = slicer;
    }

    private boolean inHmfPanel(final GeneAnnotation g) {
        return slicer.hmfRegions().stream().anyMatch(r -> r.gene().equals(g.getGeneName()));
    }

    private boolean intronicDisruption(final StructuralVariantAnnotation sv) {
        for (final GeneAnnotation g : sv.getStartAnnotations().getGenes()) {
            if (sv.getEndAnnotations()
                    .getGenes()
                    .stream()
                    .filter(o -> o.getCanonical().isIntronic() && g.getCanonical().isIntronic()
                            && o.getCanonical().getExonUpstream() == g.getCanonical().getExonUpstream())
                    .count() > 0) {
                return true;
            }
        }
        return false;
    }

    private String exonDescription(final TranscriptAnnotation t) {
        if (t.isPromoter()) {
            return "Promoter Region";
        } else if (t.isExonic()) {
            return String.format("Exon %d / %d", t.getExonUpstream(), t.getExonMax());
        } else {
            return String.format("Intron %d / %d", t.getExonUpstream(), t.getExonMax() - 1);
        }
    }

    private String exonSelection(final TranscriptAnnotation t, final boolean upstream) {
        return upstream
                ? String.format("Upstream Exon %d / %d", t.getExonUpstream(), t.getExonMax())
                : String.format("Downstream Exon %d / %d", t.getExonDownstream(), t.getExonMax());
    }

    private List<StructuralVariantAnalysis.GeneFusion> processFusions(final List<StructuralVariantAnnotation> annotations) {

        final List<Pair<TranscriptAnnotation, TranscriptAnnotation>> transcriptFusions = Lists.newArrayList();

        for (final StructuralVariantAnnotation sv : annotations) {

            final List<Pair<TranscriptAnnotation, TranscriptAnnotation>> svFusions = Lists.newArrayList();

            for (final GeneAnnotation g : sv.getStartAnnotations().getGenes()) {

                final boolean g_upstream = g.getStrand() * g.getBreakend().getOrientation() > 0;

                for (final GeneAnnotation o : sv.getEndAnnotations().getGenes()) {

                    if (g.getGeneName().equals(o.getGeneName())) {
                        continue;
                    }

                    final boolean o_upstream = o.getStrand() * o.getBreakend().getOrientation() > 0;

                    // can't both be upstream
                    if (g_upstream == o_upstream) {
                        continue;
                    }

                    for (final TranscriptAnnotation t1 : g.getTranscripts()) {
                        if (!t1.isIntronic()) {
                            continue;
                        }

                        for (final TranscriptAnnotation t2 : o.getTranscripts()) {
                            if (!t2.isIntronic()) {
                                continue;
                            }

                            if (g_upstream
                                    ? t1.getExonUpstreamPhase() == t2.getExonDownstreamPhase()
                                    : t1.getExonDownstreamPhase() == t2.getExonUpstreamPhase()) {
                                svFusions.add(Pair.of(t1, t2));
                            }

                        }
                    }
                }
            }

            Optional<Pair<TranscriptAnnotation, TranscriptAnnotation>> fusion =
                    svFusions.stream().filter(p -> p.getLeft().isCanonical() && p.getRight().isCanonical()).findFirst();

            if (fusion.isPresent()) {
                transcriptFusions.add(fusion.get());
                continue;
            }

            fusion = svFusions.stream()
                    .filter(p -> p.getLeft().isCanonical() || p.getRight().isCanonical())
                    .sorted(Comparator.comparingInt(a -> a.getLeft().getExonMax() + a.getRight().getExonMax()))
                    .reduce((a, b) -> b); // get longest

            if (fusion.isPresent()) {
                transcriptFusions.add(fusion.get());
                continue;
            }

            svFusions.stream()
                    .sorted(Comparator.comparingInt(a -> a.getLeft().getExonMax() + a.getRight().getExonMax()))
                    .reduce((a, b) -> b) // get longest
                    .ifPresent(transcriptFusions::add);
        }

        // transform results to reported details

        final List<StructuralVariantAnalysis.GeneFusion> result = Lists.newArrayList();
        for (final Pair<TranscriptAnnotation, TranscriptAnnotation> fusion : transcriptFusions) {

            final GeneAnnotation left = fusion.getLeft().getGene();
            final GeneAnnotation right = fusion.getRight().getGene();

            final StructuralVariantAnalysis.GeneFusion details = new StructuralVariantAnalysis.GeneFusion();
            details.Type = left.getBreakend().getStructuralVariant().getVariant().type().toString();

            details.GeneStart = left.getGeneName();
            details.Start = left.getBreakend().getPositionString();
            details.GeneContextStart = exonSelection(fusion.getLeft(), left.getStrand() * left.getBreakend().getOrientation() > 0);

            details.GeneEnd = right.getGeneName();
            details.End = right.getBreakend().getPositionString();
            details.GeneContextEnd = exonSelection(fusion.getRight(), right.getStrand() * right.getBreakend().getOrientation() > 0);

            result.add(details);
            annotations.remove(fusion.getLeft().getGene().getBreakend().getStructuralVariant());
        }

        return result;
    }

    private List<StructuralVariantAnalysis.GeneDisruption> processDisruptions(final List<StructuralVariantAnnotation> annotations) {

        final List<GeneAnnotation> genes = Lists.newArrayList();
        for (final StructuralVariantAnnotation sv : annotations) {

            if (intronicDisruption(sv)) {
                continue;
            }

            genes.addAll(sv.getStartAnnotations().getGenes());
            genes.addAll(sv.getEndAnnotations().getGenes());
        }

        final ArrayListMultimap<String, GeneAnnotation> geneMap = ArrayListMultimap.create();
        for (final GeneAnnotation g : genes) {
            if (!inHmfPanel(g)) {
                continue;
            }
            geneMap.put(g.getGeneName(), g);
        }

        final List<StructuralVariantAnalysis.GeneDisruption> disruptions = Lists.newArrayList();
        for (final String geneName : geneMap.keySet()) {
            for (final GeneAnnotation g : geneMap.get(geneName)) {

                final StructuralVariantAnalysis.GeneDisruption disruption = new StructuralVariantAnalysis.GeneDisruption();
                disruption.GeneName = geneName;
                disruption.Location = g.getBreakend().getPositionString();
                disruption.GeneContext = exonDescription(g.getCanonical());
                disruption.Orientation = g.getBreakend().getOrientation() > 0 ? "5\"" : "3\"";
                disruption.Partner = g.getOtherBreakend().getPositionString();
                disruption.HGVS = "TODO";
                disruption.Type = g.getBreakend().getStructuralVariant().getVariant().type().toString();
                disruption.VAF = PatientReportFormat.formatNullablePercent(g.getBreakend().getAlleleFrequency());

                disruptions.add(disruption);
                annotations.remove(g.getBreakend().getStructuralVariant());
            }
        }

        return disruptions;
    }

    @NotNull
    public StructuralVariantAnalysis run(@NotNull final List<StructuralVariant> variants) {

        final List<StructuralVariantAnnotation> annotations = annotator.annotateVariants(variants);
        final List<StructuralVariantAnalysis.GeneFusion> fusions = processFusions(annotations);
        final List<StructuralVariantAnalysis.GeneDisruption> disruptions = processDisruptions(annotations);

        return new StructuralVariantAnalysis(annotations, fusions, disruptions);
    }
}
