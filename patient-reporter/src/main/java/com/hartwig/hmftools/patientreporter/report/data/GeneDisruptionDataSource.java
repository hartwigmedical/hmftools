package com.hartwig.hmftools.patientreporter.report.data;

import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.exonDescription;
import static com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat.ploidyToCopiesString;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneDisruptionDataSource {

    private static final Logger LOGGER = LogManager.getLogger(GeneDisruptionDataSource.class);

    public static final FieldBuilder<?> CHROMOSOME_FIELD = field("chromosome", String.class);
    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> RANGE_FIELD = field("range", String.class);
    public static final FieldBuilder<?> TYPE_FIELD = field("type", String.class);
    public static final FieldBuilder<?> COPIES_FIELD = field("copies", String.class);
    public static final FieldBuilder<?> GENE_MIN_COPIES = field("gene_min_copies", String.class);
    public static final FieldBuilder<?> GENE_MAX_COPIES = field("gene_max_copies", String.class);

    private GeneDisruptionDataSource() {
    }

    @NotNull
    public static FieldBuilder<?>[] geneDisruptionFields() {
        return new FieldBuilder<?>[] { CHROMOSOME_FIELD, GENE_FIELD, RANGE_FIELD, TYPE_FIELD, COPIES_FIELD, GENE_MIN_COPIES,
                GENE_MAX_COPIES };
    }

    @NotNull
    public static JRDataSource fromGeneDisruptions(@NotNull FittedPurityStatus fitStatus, @NotNull List<GeneDisruption> disruptions) {
        final DRDataSource dataSource = new DRDataSource(CHROMOSOME_FIELD.getName(),
                GENE_FIELD.getName(),
                RANGE_FIELD.getName(),
                TYPE_FIELD.getName(),
                COPIES_FIELD.getName(),
                GENE_MIN_COPIES.getName(),
                GENE_MAX_COPIES.getName());

        Map<StructuralVariant, Pair<GeneDisruption, GeneDisruption>> pairedMap = mapDisruptionsPerStructuralVariant(disruptions);

        List<GeneDisruptionData> disruptionDataList = Lists.newArrayList();
        for (Map.Entry<StructuralVariant, Pair<GeneDisruption, GeneDisruption>> entry : pairedMap.entrySet()) {
            Pair<GeneDisruption, GeneDisruption> pairedDisruption = entry.getValue();
            GeneDisruption primaryDisruption = pairedDisruption.getLeft();

            disruptionDataList.add(ImmutableGeneDisruptionData.builder()
                    .chromosome(chromosomeField(primaryDisruption))
                    .gene(gene(primaryDisruption).geneName())
                    .type(typeField(primaryDisruption))
                    .range(rangeField(pairedDisruption))
                    .copies(copiesField(primaryDisruption))
                    .geneMinCopies(0)
                    .geneMaxCopies(0)
                    .exonUpstream(primaryDisruption.linkedAnnotation().exonUpstream())
                    .build());
        }

        for (GeneDisruptionData disruption : sort(disruptionDataList)) {
            dataSource.add(disruption.chromosome(),
                    disruption.gene(),
                    disruption.range(),
                    disruption.type(),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, disruption.copies()),
                    Strings.EMPTY,
                    Strings.EMPTY);
            // TODO (KODU): Propagate copy number info here - DEV-548
                    //                    PatientReportFormat.correctValueForFitStatus(fitStatus, String.valueOf(disruption.geneMinCopies())),
                    //                    PatientReportFormat.correctValueForFitStatus(fitStatus, String.valueOf(disruption.geneMaxCopies())));
        }

        return dataSource;
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
                left = isUpstream(disruptions.get(0)) ? disruptions.get(0) : disruptions.get(1);
                right = !isUpstream(disruptions.get(0)) ? disruptions.get(0) : disruptions.get(1);
            }

            pairedMap.put(entry.getKey(), Pair.of(left, right));
        }

        return pairedMap;
    }

    @NotNull
    private static List<GeneDisruptionData> sort(@NotNull List<GeneDisruptionData> disruptions) {
        return disruptions.stream().sorted((disruption1, disruption2) -> {
            String locationAndGene1 = zeroPrefixed(disruption1.chromosome()) + disruption1.gene();
            String locationAndGene2 = zeroPrefixed(disruption2.chromosome()) + disruption2.gene();

            if (locationAndGene1.equals(locationAndGene2)) {
                return disruption1.exonUpstream() - disruption2.exonUpstream();
            } else {
                return locationAndGene1.compareTo(locationAndGene2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String zeroPrefixed(@NotNull String chromosome) {
        // KODU: First remove q or p arm if present.
        int armStart = chromosome.indexOf("q");
        if (armStart < 0) {
            armStart = chromosome.indexOf("p");
        }

        String chromosomeForFilter = armStart > 0 ? chromosome.substring(0, armStart) : chromosome;

        try {
            int chromosomeIndex = Integer.valueOf(chromosomeForFilter);
            if (chromosomeIndex < 10) {
                return "0" + chromosome;
            } else {
                return chromosome;
            }
        } catch (NumberFormatException exception) {
            return chromosome;
        }
    }

    @NotNull
    private static String chromosomeField(@NotNull GeneDisruption disruption) {
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
    private static String typeField(@NotNull GeneDisruption disruption) {
        return gene(disruption).variant().type().name();
    }

    @NotNull
    private static String copiesField(@NotNull GeneDisruption disruption) {
        return ploidyToCopiesString(gene(disruption).variant().ploidy());
    }

    @NotNull
    private static GeneAnnotation gene(@NotNull GeneDisruption disruption) {
        return disruption.linkedAnnotation().parent();
    }

    private static boolean isUpstream(@NotNull GeneDisruption disruption) {
        GeneAnnotation gene = gene(disruption);
        return gene.variant().orientation(gene.isStart()) > 0;
    }
}
