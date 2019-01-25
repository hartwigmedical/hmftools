package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.patientreporter.loadStructuralVariants.Disruption;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableGeneDisruptionFactory {

    private static final Logger LOGGER = LogManager.getLogger(ReportableGeneDisruptionFactory.class);

    private ReportableGeneDisruptionFactory() {
    }


    @NotNull
    public static List<ReportableGeneDisruption> disruptionConvertGeneDisruption(@Nullable  List<Disruption> disruptions,
            @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        List<ReportableGeneDisruption> reportableDisruptions = Lists.newArrayList();
        Map<String, GeneCopyNumber> copyNumberPerGene = toGeneMap(geneCopyNumbers);
        if (disruptions != null) {
            for (Disruption disruption: disruptions) {
                GeneCopyNumber copyNumber = copyNumberPerGene.get(disruption.gene());
                reportableDisruptions.add(ImmutableReportableGeneDisruption.builder()
                        .location(disruption.chromosome()) // add band
                        .gene(disruption.gene())
                        .type(disruption.type())
                        .range(disruption.regionType())
                        .ploidy(0)
                        .geneMinCopies((int) Math.max(0, Math.round(copyNumber.minCopyNumber())))
                        .geneMaxCopies((int) Math.max(0, Math.round(copyNumber.maxCopyNumber())))
                     //   .firstAffectedExon(Integer.parseInt(disruption.exon()))
                        .build());
        }

        }
        return reportableDisruptions;
    }

//    @NotNull
//    public static List<ReportableGeneDisruption> toReportableGeneDisruptions(@NotNull List<GeneDisruption> disruptions,
//            @NotNull List<GeneCopyNumber> geneCopyNumbers) {
//        LOGGER.debug("Generating reportable disruptions based on {} disruptions", disruptions.size());
//        Map<String, GeneCopyNumber> copyNumberPerGene = toGeneMap(geneCopyNumbers);
//        Map<SvAndGeneKey, Pair<GeneDisruption, GeneDisruption>> pairedMap = mapDisruptionsPerStructuralVariant(disruptions);
//
//        List<ReportableGeneDisruption> reportableDisruptions = Lists.newArrayList();
//        for (Pair<GeneDisruption, GeneDisruption> pairedDisruption : pairedMap.values()) {
//            GeneDisruption primaryDisruption = pairedDisruption.getLeft();
//
//            String gene = gene(primaryDisruption).GeneName;
//            GeneCopyNumber copyNumber = copyNumberPerGene.get(gene);
//            Integer geneMinCopies = null;
//            Integer geneMaxCopies = null;
//            if (copyNumber != null) {
//                geneMinCopies = (int) Math.max(0, Math.round(copyNumber.minCopyNumber()));
//                geneMaxCopies = (int) Math.max(0, Math.round(copyNumber.maxCopyNumber()));
//            } else {
//                LOGGER.warn("Could not find copy number for disruption annotation for gene  " + gene);
//            }
//
//            Double ploidy = gene(primaryDisruption).variant().ploidy();
//
//            reportableDisruptions.add(ImmutableReportableGeneDisruption.builder()
//                    .location(locationField(primaryDisruption))
//                    .gene(gene)
//                    .type(gene(primaryDisruption).variant().type())
//                    .range(rangeField(pairedDisruption))
//                    // Not sure when ploidy ever would be null
//                    .ploidy(ploidy != null ? ploidy : Double.NaN)
//                    .geneMinCopies(geneMinCopies)
//                    .geneMaxCopies(geneMaxCopies)
//                    .firstAffectedExon(primaryDisruption.linkedAnnotation().exonUpstream())
//                    .build());
//        }
//
//        return reportableDisruptions;
//    }

    @NotNull
    private static Map<String, GeneCopyNumber> toGeneMap(@NotNull List<GeneCopyNumber> geneCopyNumbers) {
        Map<String, GeneCopyNumber> copyNumberPerGeneMap = Maps.newHashMap();
        for (GeneCopyNumber copyNumber : geneCopyNumbers) {
            copyNumberPerGeneMap.put(copyNumber.gene(), copyNumber);
        }
        return copyNumberPerGeneMap;
    }

    @NotNull
    private static Map<SvAndGeneKey, Pair<GeneDisruption, GeneDisruption>> mapDisruptionsPerStructuralVariant(
            @NotNull List<GeneDisruption> disruptions) {
        Map<SvAndGeneKey, List<GeneDisruption>> disruptionsPerSvAndGene = Maps.newHashMap();
        for (GeneDisruption disruption : disruptions) {
            StructuralVariant variant = disruption.linkedAnnotation().parent().variant();
            String gene = disruption.linkedAnnotation().parent().GeneName;
            SvAndGeneKey key = new SvAndGeneKey(variant, gene);
            List<GeneDisruption> currentDisruptions = disruptionsPerSvAndGene.get(key);
            if (currentDisruptions == null) {
                currentDisruptions = Lists.newArrayList();
            }
            currentDisruptions.add(disruption);
            disruptionsPerSvAndGene.put(key, currentDisruptions);
        }

        return toPairedMap(disruptionsPerSvAndGene);
    }

    @NotNull
    private static Map<SvAndGeneKey, Pair<GeneDisruption, GeneDisruption>> toPairedMap(
            @NotNull Map<SvAndGeneKey, List<GeneDisruption>> disruptionsPerVariant) {
        Map<SvAndGeneKey, Pair<GeneDisruption, GeneDisruption>> pairedMap = Maps.newHashMap();

        for (Map.Entry<SvAndGeneKey, List<GeneDisruption>> entry : disruptionsPerVariant.entrySet()) {
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
                boolean firstBeforeSecond =
                        disruptions.get(0).linkedAnnotation().exonUpstream() <= disruptions.get(1).linkedAnnotation().exonUpstream();
                left = firstBeforeSecond ? disruptions.get(0) : disruptions.get(1);
                right = firstBeforeSecond ? disruptions.get(1) : disruptions.get(0);
            }

            pairedMap.put(entry.getKey(), Pair.of(left, right));
        }

        return pairedMap;
    }

    @NotNull
    private static String locationField(@NotNull GeneDisruption disruption) {
        GeneAnnotation gene = gene(disruption);
        return gene.chromosome() + gene.karyotypeBand();
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
        GeneAnnotation gene = gene(disruption);
        return gene.orientation() * gene.Strand < 0;
    }

    private static class SvAndGeneKey {
        @NotNull
        private final StructuralVariant variant;
        @NotNull
        private final String gene;

        private SvAndGeneKey(@NotNull final StructuralVariant variant, @NotNull final String gene) {
            this.variant = variant;
            this.gene = gene;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final SvAndGeneKey that = (SvAndGeneKey) o;
            return Objects.equals(variant, that.variant) && Objects.equals(gene, that.gene);
        }

        @Override
        public int hashCode() {
            return Objects.hash(variant, gene);
        }
    }
}
