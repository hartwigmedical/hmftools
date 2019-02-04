package com.hartwig.hmftools.patientreporter.structural;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ReportableGeneDisruptionFactory {

    private static final Logger LOGGER = LogManager.getLogger(ReportableGeneDisruptionFactory.class);

    private ReportableGeneDisruptionFactory() {
    }

    @NotNull
    public static List<ReportableGeneDisruption> disruptionConvertGeneDisruption(@NotNull List<Disruption> disruptions,
            @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        LOGGER.debug("Generating reportable disruptions based on {} disruptions", disruptions.size());

        List<ReportableGeneDisruption> reportableDisruptions = Lists.newArrayList();
        Map<String, GeneCopyNumber> copyNumberPerGene = toGeneMap(geneCopyNumbers);
        Map<SvAndGeneKey, Pair<Disruption, Disruption>> pairedMap = mapDisruptionsPerStructuralVariant(disruptions);

        for (Pair<Disruption, Disruption> pairedDisruption : pairedMap.values()) {
            Disruption primaryDisruption = pairedDisruption.getLeft();

            GeneCopyNumber copyNumber = copyNumberPerGene.get(primaryDisruption.gene());
            reportableDisruptions.add(ImmutableReportableGeneDisruption.builder()
                    .location(primaryDisruption.chromosome() + primaryDisruption.chrBand())
                    .gene(primaryDisruption.gene())
                    .type(primaryDisruption.type())
                    .range(rangeField(pairedDisruption))
                    .ploidy(primaryDisruption.ploidy())
                    .geneMinCopies((int) Math.max(0, Math.round(copyNumber.minCopyNumber())))
                    .geneMaxCopies((int) Math.max(0, Math.round(copyNumber.maxCopyNumber())))
                    .firstAffectedExon(primaryDisruption.exonUp())
                    .build());
        }

        return reportableDisruptions;
    }

    @NotNull
    private static Map<String, GeneCopyNumber> toGeneMap(@NotNull List<GeneCopyNumber> geneCopyNumbers) {
        Map<String, GeneCopyNumber> copyNumberPerGeneMap = Maps.newHashMap();
        for (GeneCopyNumber copyNumber : geneCopyNumbers) {
            copyNumberPerGeneMap.put(copyNumber.gene(), copyNumber);
        }
        return copyNumberPerGeneMap;
    }

    @NotNull
    private static Map<SvAndGeneKey, Pair<Disruption, Disruption>> mapDisruptionsPerStructuralVariant(
            @NotNull List<Disruption> disruptions) {
        Map<SvAndGeneKey, List<Disruption>> disruptionsPerSvAndGene = Maps.newHashMap();
        for (Disruption disruption : disruptions) {
            SvAndGeneKey key = new SvAndGeneKey(disruption.svId(), disruption.gene());
            List<Disruption> currentDisruptions = disruptionsPerSvAndGene.get(key);
            if (currentDisruptions == null) {
                currentDisruptions = Lists.newArrayList();
            }
            currentDisruptions.add(disruption);
            disruptionsPerSvAndGene.put(key, currentDisruptions);
        }

        return toPairedMap(disruptionsPerSvAndGene);
    }

    @NotNull
    private static Map<SvAndGeneKey, Pair<Disruption, Disruption>> toPairedMap(
            @NotNull Map<SvAndGeneKey, List<Disruption>> disruptionsPerVariant) {
        Map<SvAndGeneKey, Pair<Disruption, Disruption>> pairedMap = Maps.newHashMap();

        for (Map.Entry<SvAndGeneKey, List<Disruption>> entry : disruptionsPerVariant.entrySet()) {
            List<Disruption> disruptions = entry.getValue();

            if (disruptions.size() != 1 && disruptions.size() != 2) {
                LOGGER.warn("Found unusual number of disruptions on single event: " + disruptions.size());
                continue;
            }

            Disruption left;
            Disruption right;
            if (disruptions.size() == 1) {
                left = disruptions.get(0);
                right = null;
            } else {
                boolean firstBeforeSecond = disruptions.get(0).exonUp() <= disruptions.get(1).exonUp();
                left = firstBeforeSecond ? disruptions.get(0) : disruptions.get(1);
                right = firstBeforeSecond ? disruptions.get(1) : disruptions.get(0);
            }

            pairedMap.put(entry.getKey(), Pair.of(left, right));
        }

        return pairedMap;
    }

    @NotNull
    private static String rangeField(@NotNull Pair<Disruption, Disruption> pairedDisruption) {
        Disruption primary = pairedDisruption.getLeft();
        Disruption secondary = pairedDisruption.getRight();

        if (secondary == null) {
            return exonDescription(primary.exonUp(), primary.exonDown()) + (isUpstream(primary) ? " Upstream" : " Downstream");
        } else {
            return exonDescription(primary.exonUp(), primary.exonDown()) + " -> " + exonDescription(secondary.exonUp(),
                    secondary.exonDown());
        }
    }

    private static boolean isUpstream(@NotNull Disruption disruption) {
        return disruption.orientation() * disruption.strand() < 0;
    }

    private static class SvAndGeneKey {
        @NotNull
        private final String variantId;
        @NotNull
        private final String gene;

        private SvAndGeneKey(@NotNull final String variantId, @NotNull final String gene) {
            this.variantId = variantId;
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
            return Objects.equals(variantId, that.variantId) && Objects.equals(gene, that.gene);
        }

        @Override
        public int hashCode() {

            return Objects.hash(variantId, gene);
        }
    }
}
