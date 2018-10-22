package com.hartwig.hmftools.patientreporter.disruption;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ReportableGeneDisruptionFactory {

    private static final Logger LOGGER = LogManager.getLogger(ReportableGeneDisruptionFactory.class);

    private ReportableGeneDisruptionFactory() {
    }

    @NotNull
    public static List<ReportableGeneDisruption> toReportableGeneDisruptions(@NotNull List<GeneDisruption> disruptions) {
        Map<StructuralVariant, Pair<GeneDisruption, GeneDisruption>> pairedMap = mapDisruptionsPerStructuralVariant(disruptions);

        List<ReportableGeneDisruption> reportableDisruptions = Lists.newArrayList();
        for (Map.Entry<StructuralVariant, Pair<GeneDisruption, GeneDisruption>> entry : pairedMap.entrySet()) {
            Pair<GeneDisruption, GeneDisruption> pairedDisruption = entry.getValue();
            GeneDisruption primaryDisruption = pairedDisruption.getLeft();

            Double ploidy = gene(primaryDisruption).variant().ploidy();

            reportableDisruptions.add(ImmutableReportableGeneDisruption.builder()
                    .location(locationField(primaryDisruption))
                    .gene(gene(primaryDisruption).geneName())
                    .type(gene(primaryDisruption).variant().type())
                    .range(rangeField(pairedDisruption))
                    // KODU: Not sure when ploidy ever would be null
                    .ploidy(ploidy != null ? ploidy : Double.NaN)
                    .geneMinCopies(0)
                    .geneMaxCopies(0)
                    .firstAffectedExon(primaryDisruption.linkedAnnotation().exonUpstream())
                    .build());
        }

        return reportableDisruptions;
    }

    @NotNull
    private static Map<StructuralVariant, Pair<GeneDisruption, GeneDisruption>> mapDisruptionsPerStructuralVariant(
            @NotNull List<GeneDisruption> disruptions) {
        Map<StructuralVariant, List<GeneDisruption>> disruptionsPerVariant = Maps.newHashMap();
        for (GeneDisruption disruption : disruptions) {
            StructuralVariant variant = disruption.linkedAnnotation().parent().variant();
            List<GeneDisruption> currentDisruptions = disruptionsPerVariant.get(variant);
            if (currentDisruptions == null) {
                currentDisruptions = Lists.newArrayList();
            }
            currentDisruptions.add(disruption);
            disruptionsPerVariant.put(variant, currentDisruptions);
        }

        return toPairedMap(disruptionsPerVariant);
    }

    @NotNull
    private static Map<StructuralVariant, Pair<GeneDisruption, GeneDisruption>> toPairedMap(
            @NotNull Map<StructuralVariant, List<GeneDisruption>> disruptionsPerVariant) {
        Map<StructuralVariant, Pair<GeneDisruption, GeneDisruption>> pairedMap = Maps.newHashMap();

        for (Map.Entry<StructuralVariant, List<GeneDisruption>> entry : disruptionsPerVariant.entrySet()) {
            List<GeneDisruption> disruptions = entry.getValue();

            if (disruptions.size() != 1 && disruptions.size() != 2) {
                LOGGER.warn("Found unusual number of disruptions on single event: " + disruptions.size());
                continue;
            }

            GeneDisruption left;
            GeneDisruption right;
            if (disruptions.size() == 1) {
                left = disruptions.get(0);
                right = null;
            } else {
                // TODO (KODU): Should I use isStart/isEnd rather than upstream/downstream?
                left = isUpstream(disruptions.get(0)) ? disruptions.get(0) : disruptions.get(1);
                right = !isUpstream(disruptions.get(0)) ? disruptions.get(0) : disruptions.get(1);
            }

            pairedMap.put(entry.getKey(), Pair.of(left, right));
        }

        return pairedMap;
    }

    @NotNull
    private static String locationField(@NotNull GeneDisruption disruption) {
        GeneAnnotation gene = gene(disruption);
        return gene.variant().chromosome(gene.isStart()) + gene.karyotypeBand();
    }

    @NotNull
    private static String rangeField(@NotNull Pair<GeneDisruption, GeneDisruption> pairedDisruption) {
        GeneDisruption primary = pairedDisruption.getLeft();
        GeneDisruption secondary = pairedDisruption.getRight();
        if (secondary == null) {
            return exonDescription(primary.linkedAnnotation()) + (isUpstream(primary) ? " Upstream" : " Downstream");
        } else {
            return exonDescription(primary.linkedAnnotation()) + " -> " + exonDescription(secondary.linkedAnnotation());
        }
    }

    @NotNull
    private static GeneAnnotation gene(@NotNull GeneDisruption disruption) {
        return disruption.linkedAnnotation().parent();
    }

    private static boolean isUpstream(@NotNull GeneDisruption disruption) {
        // KODU (TODO): Figure out whether definition is correct!
        GeneAnnotation gene = gene(disruption);
        return gene.variant().orientation(gene.isStart()) > 0;
    }
}
