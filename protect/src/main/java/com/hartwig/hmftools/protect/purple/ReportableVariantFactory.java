package com.hartwig.hmftools.protect.purple;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableVariantFactory {

    private ReportableVariantFactory() {
    }

    @NotNull
    public static List<ReportableVariant> reportableGermlineVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> germlineDriverCatalog, boolean hasReliablePurity) {
        return reportableVariants(variants, germlineDriverCatalog, hasReliablePurity, ReportableVariantSource.GERMLINE);
    }

    @NotNull
    public static List<ReportableVariant> reportableSomaticVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalog, boolean hasReliablePurity) {
        List<DriverCatalog> mutationCatalog =
                driverCatalog.stream().filter(x -> x.driver() == DriverType.MUTATION).collect(Collectors.toList());
        return reportableVariants(variants, mutationCatalog, hasReliablePurity, ReportableVariantSource.SOMATIC);

    }

    @NotNull
    private static List<ReportableVariant> reportableVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalog, boolean hasReliablePurity, ReportableVariantSource source) {
        Map<String, DriverCatalog> mutationDriverMap = driverCatalog.stream().collect(Collectors.toMap(DriverCatalog::gene, x -> x));
        Map<String, HmfTranscriptRegion> genePanel = HmfGenePanelSupplier.allGenesMap37();

        List<ReportableVariant> result = Lists.newArrayList();
        for (SomaticVariant variant : variants) {
            if (variant.reported()) {
                DriverCatalog geneDriver = mutationDriverMap.get(variant.gene());
                assert geneDriver != null;

                HmfTranscriptRegion canonicalTranscript = genePanel.get(variant.gene());

                ReportableVariant reportable =
                        fromVariant(variant, hasReliablePurity, source, canonicalTranscript).driverLikelihood(geneDriver.driverLikelihood()).build();
                result.add(reportable);
            }
        }

        return result;
    }

    @NotNull
    public static List<ReportableVariant> mergeVariantLists(@NotNull List<ReportableVariant> list1,
            @NotNull List<ReportableVariant> list2) {
        List<ReportableVariant> result = Lists.newArrayList();

        Map<String, Double> maxLikelihood = Maps.newHashMap();
        for (ReportableVariant variant : list1) {
            maxLikelihood.merge(variant.gene(), variant.driverLikelihood(), Math::max);
        }

        for (ReportableVariant variant : list2) {
            maxLikelihood.merge(variant.gene(), variant.driverLikelihood(), Math::max);
        }

        for (ReportableVariant variant : list1) {
            result.add(ImmutableReportableVariant.builder().from(variant).driverLikelihood(maxLikelihood.get(variant.gene())).build());
        }

        for (ReportableVariant variant : list2) {
            result.add(ImmutableReportableVariant.builder().from(variant).driverLikelihood(maxLikelihood.get(variant.gene())).build());
        }

        return result;
    }

    @NotNull
    private static ImmutableReportableVariant.Builder fromVariant(@NotNull SomaticVariant variant, boolean hasReliablePurity,
            ReportableVariantSource source, @Nullable HmfTranscriptRegion canonicalTranscript) {
        return ImmutableReportableVariant.builder()
                .type(variant.type())
                .source(source)
                .gene(variant.gene())
                .chromosome(variant.chromosome())
                .position(variant.position())
                .ref(variant.ref())
                .alt(variant.alt())
                .canonicalTranscript(canonicalTranscript == null ? Strings.EMPTY : canonicalTranscript.geneID())
                .canonicalCodingEffect(variant.canonicalCodingEffect())
                .canonicalHgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .canonicalHgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .totalCopyNumber(variant.adjustedCopyNumber())
                .copyNumber(copyNumberString(variant.adjustedCopyNumber(), hasReliablePurity))
                .tVafString(vafString(calcAlleleCopyNumber(variant.adjustedCopyNumber(), variant.adjustedVAF()),
                        variant.adjustedCopyNumber(),
                        hasReliablePurity))
                .alleleCopyNumber(calcAlleleCopyNumber(variant.adjustedCopyNumber(), variant.adjustedVAF()))
                .hotspot(variant.hotspot())
                .clonalLikelihood(variant.clonalLikelihood())
                .biallelic(variant.biallelic());
    }

    private static double calcAlleleCopyNumber(double adjustedCopyNumber, double adjustedVAF) {
        return adjustedCopyNumber * Math.max(0, Math.min(1, adjustedVAF));
    }

    @NotNull
    public static String copyNumberString(double copyNumber, boolean hasReliablePurity) {
        if (!hasReliablePurity) {
            return DataUtil.NA_STRING;
        }

        return String.valueOf(Math.round(Math.max(0, copyNumber)));
    }

    @NotNull
    public static String vafString(double alleleCopyNumber, double totalCopyNumber, boolean hasReliablePurity) {
        if (!hasReliablePurity) {
            return DataUtil.NA_STRING;
        }
        double vaf = alleleCopyNumber / totalCopyNumber;

        return DataUtil.formatPercentage(100 * Math.max(0, Math.min(1, vaf)));
    }
}
