package com.hartwig.hmftools.patientreporter.variants;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.slicing.HmfSlicer;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;
import com.hartwig.hmftools.svannotation.GeneAnnotation;
import com.hartwig.hmftools.svannotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.svannotation.StructuralVariantAnnotator;
import com.hartwig.hmftools.svannotation.TranscriptAnnotation;

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
            return String.format("In exon %d / %d", t.getExonUpstream(), t.getExonMax());
        } else {
            return String.format("Between exon %d and %d / %d", t.getExonUpstream(), t.getExonDownstream(), t.getExonMax());
        }
    }

    private List<StructuralVariantAnalysis.GeneFusion> processFusions(final List<StructuralVariantAnnotation> annotations) {
        return Collections.emptyList();
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
                disruption.Partner = g.getOtherBreakend().getPositionString();
                disruption.HGVS = "TODO";
                disruption.Type = g.getBreakend().getStructuralVariant().getVariant().type().toString();
                disruption.Orientation = g.getBreakend().getOrientation() > 0 ? "5\"" : "3\"";
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
