package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogKey;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogMap;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableVariantFactory
{
    private static final Logger LOGGER = LogManager.getLogger(ReportableVariantFactory.class);

    private ReportableVariantFactory()
    {
    }

    @NotNull
    public static List<ReportableVariant> toReportableGermlineVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> germlineDriverCatalog)
    {
        List<DriverCatalog> germlineMutationCatalog =
                germlineDriverCatalog.stream().filter(x -> x.driver() == DriverType.GERMLINE_MUTATION).collect(Collectors.toList());
        return toReportableVariants(variants, germlineMutationCatalog, ReportableVariantSource.GERMLINE);
    }

    @NotNull
    public static List<ReportableVariant> toReportableSomaticVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> somaticDriverCatalog)
    {
        List<DriverCatalog> somaticMutationCatalog =
                somaticDriverCatalog.stream().filter(x -> x.driver() == DriverType.MUTATION).collect(Collectors.toList());
        return toReportableVariants(variants, somaticMutationCatalog, ReportableVariantSource.SOMATIC);
    }

    @NotNull
    private static List<ReportableVariant> toReportableVariants(@NotNull List<SomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalog, @NotNull ReportableVariantSource source)
    {
        Map<DriverCatalogKey, DriverCatalog> driverMap = DriverCatalogMap.toDriverMap(driverCatalog);
        List<ReportableVariant> result = Lists.newArrayList();

        for(SomaticVariant variant : variants)
        {
            if(variant.reported())
            {
                ImmutableReportableVariant.Builder builder = fromVariant(variant, source);

                DriverCatalog canonicalDriver = findCanonicalEntryForVariant(driverMap, variant);
                if(canonicalDriver != null)
                {
                    result.add(builder.driverLikelihood(canonicalDriver.driverLikelihood())
                            .transcript(canonicalDriver.transcript())
                            .isCanonical(canonicalDriver.isCanonical())
                            .build());
                }

                DriverCatalog nonCanonicalDriver = findNonCanonicalEntryForVariant(driverMap, variant);
                if(nonCanonicalDriver != null)
                {
                    result.add(builder.driverLikelihood(nonCanonicalDriver.driverLikelihood())
                            .transcript(nonCanonicalDriver.transcript())
                            .isCanonical(nonCanonicalDriver.isCanonical())
                            .canonicalHgvsProteinImpact(AltTranscriptReportableInfo.firstOtherHgvsProteinImpact(variant.otherReportedEffects()))
                            .canonicalHgvsCodingImpact(AltTranscriptReportableInfo.firstOtherHgvsCodingImpact(variant.otherReportedEffects()))
                            .build());
                }
            }

        }
        return result;
    }

    @Nullable
    private static DriverCatalog findCanonicalEntryForVariant(@NotNull Map<DriverCatalogKey, DriverCatalog> entries,
            @NotNull SomaticVariant variant)
    {
        assert variant.reported();

        for(DriverCatalog catalog : entries.values())
        {
            if(variant.gene().equals(catalog.gene()) && catalog.isCanonical())
            {
                if(variant.canonicalTranscript().equals(catalog.transcript()))
                {
                    return entries.get(DriverCatalogKey.create(variant.gene(), variant.canonicalTranscript()));
                }
                else
                {
                    LOGGER.warn("Canonical driver entry on transcript {} does not match canonical transcript {} on variant",
                            catalog.transcript(),
                            variant.canonicalTranscript());
                }
            }
        }

        LOGGER.warn("No canonical entry found in driver catalog for gene {}", variant.gene());
        return null;
    }

    @Nullable
    private static DriverCatalog findNonCanonicalEntryForVariant(@NotNull Map<DriverCatalogKey, DriverCatalog> entries,
            @NotNull SomaticVariant variant)
    {
        assert variant.reported();

        String nonCanonicalTranscript = AltTranscriptReportableInfo.firstOtherTranscript(variant.otherReportedEffects());

        for(DriverCatalog catalog : entries.values())
        {
            if(variant.gene().equals(catalog.gene()) && !catalog.isCanonical())
            {
                if(nonCanonicalTranscript.equals(catalog.transcript()))
                {
                    return entries.get(DriverCatalogKey.create(variant.gene(), nonCanonicalTranscript));
                }
                else
                {
                    LOGGER.warn("Unexpected transcript {} for gene {} in driver catalog", catalog.transcript(), catalog.gene());
                }
            }
        }

        return null;
    }

    @NotNull
    public static List<ReportableVariant> mergeVariantLists(@NotNull List<ReportableVariant> list1,
            @NotNull List<ReportableVariant> list2)
    {
        List<ReportableVariant> result = Lists.newArrayList();

        Map<String, Double> maxLikelihoodPerGene = Maps.newHashMap();
        for(ReportableVariant variant : list1)
        {
            maxLikelihoodPerGene.merge(variant.gene(), variant.driverLikelihood(), Math::max);
        }

        for(ReportableVariant variant : list2)
        {
            maxLikelihoodPerGene.merge(variant.gene(), variant.driverLikelihood(), Math::max);
        }

        for(ReportableVariant variant : list1)
        {
            result.add(ImmutableReportableVariant.builder()
                    .from(variant)
                    .driverLikelihood(maxLikelihoodPerGene.get(variant.gene()))
                    .build());
        }

        for(ReportableVariant variant : list2)
        {
            result.add(ImmutableReportableVariant.builder()
                    .from(variant)
                    .driverLikelihood(maxLikelihoodPerGene.get(variant.gene()))
                    .build());
        }

        return result;
    }

    @NotNull
    public static ImmutableReportableVariant.Builder fromVariant(@NotNull SomaticVariant variant, @NotNull ReportableVariantSource source)
    {
        AllelicDepth rnaDepth = variant.rnaDepth();
        return ImmutableReportableVariant.builder()
                .type(variant.type())
                .source(source)
                .gene(variant.gene())
                .chromosome(variant.chromosome())
                .position(variant.position())
                .ref(variant.ref())
                .alt(variant.alt())
                .otherReportedEffects(variant.otherReportedEffects())
                .canonicalTranscript(variant.canonicalTranscript())
                .canonicalEffect(variant.canonicalEffect())
                .canonicalCodingEffect(variant.canonicalCodingEffect())
                .canonicalHgvsCodingImpact(variant.canonicalHgvsCodingImpact())
                .canonicalHgvsProteinImpact(variant.canonicalHgvsProteinImpact())
                .totalReadCount(variant.totalReadCount())
                .alleleReadCount(variant.alleleReadCount())
                .rnaAlleleReadCount(rnaDepth != null ? rnaDepth.alleleReadCount() : null)
                .rnaTotalReadCount(rnaDepth != null ? rnaDepth.totalReadCount() : null)
                .totalCopyNumber(variant.adjustedCopyNumber())
                .minorAlleleCopyNumber(variant.minorAlleleCopyNumber())
                .alleleCopyNumber(calcAlleleCopyNumber(variant.adjustedCopyNumber(), variant.adjustedVAF()))
                .hotspot(variant.hotspot())
                .clonalLikelihood(variant.clonalLikelihood())
                .biallelic(variant.biallelic())
                .genotypeStatus(variant.genotypeStatus())
                .localPhaseSet(variant.topLocalPhaseSet());
    }

    private static double calcAlleleCopyNumber(double adjustedCopyNumber, double adjustedVAF)
    {
        return adjustedCopyNumber * Math.max(0, Math.min(1, adjustedVAF));
    }
}
