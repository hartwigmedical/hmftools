package com.hartwig.hmftools.patientreporter.structural;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.fusion.ReportableDisruption;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class ReportableGeneDisruptionFactory {

    private static final Logger LOGGER = LogManager.getLogger(ReportableGeneDisruptionFactory.class);

    private ReportableGeneDisruptionFactory() {
    }

    @NotNull
    static List<ReportableGeneDisruption> convert(@NotNull List<ReportableDisruption> disruptions) {
        List<ReportableGeneDisruption> reportableDisruptions = Lists.newArrayList();
        Map<SvAndGeneKey, Pair<ReportableDisruption, ReportableDisruption>> pairedMap = mapDisruptionsPerStructuralVariant(disruptions);

        for (Pair<ReportableDisruption, ReportableDisruption> pairedDisruption : pairedMap.values()) {
            ReportableDisruption primaryDisruptionLeft = pairedDisruption.getLeft();
            ReportableDisruption primaryDisruptionRight = pairedDisruption.getRight();

            if (primaryDisruptionRight != null) {
                double lowestUndisruptedCopyNumber =
                        Math.min(primaryDisruptionLeft.undisruptedCopyNumber(), primaryDisruptionRight.undisruptedCopyNumber());

                Double ploidyLeft = primaryDisruptionLeft.junctionCopyNumber();
                Double ploidyRight = primaryDisruptionRight.junctionCopyNumber();
                if (ploidyLeft != null && ploidyRight != null && !Doubles.equal(ploidyLeft, ploidyRight)) {
                    LOGGER.warn("The disrupted copy number of a paired sv is not the same on {}", primaryDisruptionLeft.gene());
                }
                reportableDisruptions.add(ImmutableReportableGeneDisruption.builder()
                        .location(primaryDisruptionLeft.chromosome() + primaryDisruptionLeft.chrBand())
                        .gene(primaryDisruptionLeft.gene())
                        .type(primaryDisruptionLeft.type())
                        .range(rangeField(pairedDisruption))
                        .ploidy(primaryDisruptionLeft.junctionCopyNumber())
                        .undisruptedCopyNumber(Math.max(0, lowestUndisruptedCopyNumber))
                        .firstAffectedExon(primaryDisruptionLeft.exonUp())
                        .build());
            } else {
                reportableDisruptions.add(ImmutableReportableGeneDisruption.builder()
                        .location(primaryDisruptionLeft.chromosome() + primaryDisruptionLeft.chrBand())
                        .gene(primaryDisruptionLeft.gene())
                        .type(primaryDisruptionLeft.type())
                        .range(rangeField(pairedDisruption))
                        .ploidy(primaryDisruptionLeft.junctionCopyNumber())
                        .undisruptedCopyNumber(Math.max(0, primaryDisruptionLeft.undisruptedCopyNumber()))
                        .firstAffectedExon(primaryDisruptionLeft.exonUp())
                        .build());
            }
        }

        LOGGER.debug("Generated {} reportable disruptions based on {} disruptions", reportableDisruptions.size(), disruptions.size());
        return reportableDisruptions;
    }

    @NotNull
    private static Map<SvAndGeneKey, Pair<ReportableDisruption, ReportableDisruption>> mapDisruptionsPerStructuralVariant(
            @NotNull List<ReportableDisruption> disruptions) {
        Map<SvAndGeneKey, List<ReportableDisruption>> disruptionsPerSvAndGene = Maps.newHashMap();
        for (ReportableDisruption disruption : disruptions) {
            SvAndGeneKey key = new SvAndGeneKey(disruption.svId(), disruption.gene());
            List<ReportableDisruption> currentDisruptions = disruptionsPerSvAndGene.get(key);
            if (currentDisruptions == null) {
                currentDisruptions = Lists.newArrayList();
            }
            currentDisruptions.add(disruption);
            disruptionsPerSvAndGene.put(key, currentDisruptions);
        }

        return toPairedMap(disruptionsPerSvAndGene);
    }

    @NotNull
    private static Map<SvAndGeneKey, Pair<ReportableDisruption, ReportableDisruption>> toPairedMap(
            @NotNull Map<SvAndGeneKey, List<ReportableDisruption>> disruptionsPerVariant) {
        Map<SvAndGeneKey, Pair<ReportableDisruption, ReportableDisruption>> pairedMap = Maps.newHashMap();

        for (Map.Entry<SvAndGeneKey, List<ReportableDisruption>> entry : disruptionsPerVariant.entrySet()) {
            List<ReportableDisruption> disruptions = entry.getValue();

            if (disruptions.size() != 1 && disruptions.size() != 2) {
                LOGGER.warn("Found unusual number of disruptions on single event: {}", disruptions.size());
                continue;
            }

            ReportableDisruption left;
            ReportableDisruption right;
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
    private static String rangeField(@NotNull Pair<ReportableDisruption, ReportableDisruption> pairedDisruption) {
        ReportableDisruption primary = pairedDisruption.getLeft();
        ReportableDisruption secondary = pairedDisruption.getRight();

        if (secondary == null) {
            return exonDescription(primary.exonUp(), primary.exonDown()) + (isUpstream(primary) ? " Upstream" : " Downstream");
        } else {
            return exonDescription(primary.exonUp(), primary.exonDown()) + " -> " + exonDescription(secondary.exonUp(),
                    secondary.exonDown());
        }
    }

    @NotNull
    private static String exonDescription(int exonUp, int exonDown) {
        if (exonUp > 0) {
            if (exonUp == exonDown) {
                return String.format("Exon %d", exonUp);
            } else if (exonDown - exonUp == 1) {
                return String.format("Intron %d", exonUp);
            }
        } else if (exonUp == 0 && (exonDown == 1 || exonDown == 2)) {
            return "Promoter Region";
        }

        return String.format("ERROR up=%d, down=%d", exonUp, exonDown);
    }

    private static boolean isUpstream(@NotNull ReportableDisruption disruption) {
        return disruption.orientation() * disruption.strand() < 0;
    }

    private static class SvAndGeneKey {
        private final int variantId;
        @NotNull
        private final String gene;

        private SvAndGeneKey(final int variantId, @NotNull final String gene) {
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
