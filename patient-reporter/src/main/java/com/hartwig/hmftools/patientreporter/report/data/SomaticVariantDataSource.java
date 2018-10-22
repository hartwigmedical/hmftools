package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.patientreporter.algo.GeneModel;
import com.hartwig.hmftools.patientreporter.report.util.PatientReportFormat;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class SomaticVariantDataSource {

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> VARIANT_FIELD = field("variant", String.class);
    public static final FieldBuilder<?> IMPACT_FIELD = field("impact", String.class);
    public static final FieldBuilder<?> READ_DEPTH_FIELD = field("read_depth", String.class);
    public static final FieldBuilder<?> IS_HOTSPOT_FIELD = field("is_hotspot", String.class);
    public static final FieldBuilder<?> PLOIDY_VAF_FIELD = field("ploidy_vaf", String.class);
    public static final FieldBuilder<?> CLONAL_STATUS_FIELD = field("clonal_status", String.class);
    public static final FieldBuilder<?> BIALLELIC_FIELD = field("biallelic", String.class);
    public static final FieldBuilder<?> DRIVER_FIELD = field("driver", String.class);

    private SomaticVariantDataSource() {
    }

    @NotNull
    public static FieldBuilder<?>[] fields() {
        return new FieldBuilder<?>[] { GENE_FIELD, VARIANT_FIELD, READ_DEPTH_FIELD, IS_HOTSPOT_FIELD, PLOIDY_VAF_FIELD, CLONAL_STATUS_FIELD,
                BIALLELIC_FIELD, DRIVER_FIELD };
    }

    @NotNull
    public static JRDataSource fromVariants(@NotNull FittedPurityStatus fitStatus, @NotNull List<EnrichedSomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalogList, @NotNull GeneModel panelGeneModel) {
        final DRDataSource variantDataSource = new DRDataSource(GENE_FIELD.getName(),
                VARIANT_FIELD.getName(),
                IMPACT_FIELD.getName(),
                READ_DEPTH_FIELD.getName(),
                IS_HOTSPOT_FIELD.getName(),
                PLOIDY_VAF_FIELD.getName(),
                CLONAL_STATUS_FIELD.getName(),
                BIALLELIC_FIELD.getName(),
                DRIVER_FIELD.getName());

        for (final EnrichedSomaticVariant variant : sort(variants, driverCatalogList)) {
            DriverCategory driverCategory = panelGeneModel.geneDriverCategory(variant.gene());

            String displayGene =
                    panelGeneModel.drupActionableGenes().keySet().contains(variant.gene()) ? variant.gene() + " *" : variant.gene();

            String biallelic = Strings.EMPTY;
            if (driverCategory != DriverCategory.ONCO) {
                biallelic = variant.biallelic() ? "Yes" : "No";
            }

            String ploidyVaf =
                    PatientReportFormat.ploidyVafField(variant.adjustedCopyNumber(), variant.minorAllelePloidy(), variant.adjustedVAF());

            variantDataSource.add(displayGene,
                    variant.canonicalHgvsCodingImpact(),
                    variant.canonicalHgvsProteinImpact(),
                    PatientReportFormat.readDepthField(variant),
                    hotspotField(variant),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, ploidyVaf),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, clonalityField(variant)),
                    PatientReportFormat.correctValueForFitStatus(fitStatus, biallelic),
                    driverField(catalogEntryForVariant(driverCatalogList, variant)));
        }

        return variantDataSource;
    }

    @NotNull
    private static List<EnrichedSomaticVariant> sort(@NotNull List<EnrichedSomaticVariant> variants,
            @NotNull List<DriverCatalog> driverCatalogList) {
        return variants.stream().sorted((variant1, variant2) -> {
            DriverCatalog entryForVariant1 = catalogEntryForVariant(driverCatalogList, variant1);
            DriverCatalog entryForVariant2 = catalogEntryForVariant(driverCatalogList, variant2);

            // KODU: Force any hgvsCodingImpact outside of driver catalog to the bottom of table.
            double driverLikelihood1 = entryForVariant1 != null ? entryForVariant1.driverLikelihood() : -1;
            double driverLikelihood2 = entryForVariant2 != null ? entryForVariant2.driverLikelihood() : -1;
            if (Math.abs(driverLikelihood1 - driverLikelihood2) > 0.001) {
                return (driverLikelihood1 - driverLikelihood2) < 0 ? 1 : -1;
            } else {
                if (variant1.gene().equals(variant2.gene())) {
                    return variant1.canonicalHgvsCodingImpact().compareTo(variant2.canonicalHgvsCodingImpact());
                } else {
                    return variant1.gene().compareTo(variant2.gene());
                }
            }
        }).collect(Collectors.toList());
    }

    @Nullable
    private static DriverCatalog catalogEntryForVariant(@NotNull List<DriverCatalog> driverCatalogList, @NotNull SomaticVariant variant) {
        for (DriverCatalog entry : driverCatalogList) {
            if (entry.gene().equals(variant.gene())) {
                return entry;
            }
        }
        return null;
    }

    @NotNull
    private static String driverField(@Nullable DriverCatalog catalogEntry) {
        if (catalogEntry == null) {
            return Strings.EMPTY;
        }

        if (catalogEntry.driverLikelihood() > 0.8) {
            return "Likely";
        } else if (catalogEntry.driverLikelihood() > 0.2) {
            return "Potentially";
        } else {
            return "Unlikely";
        }
    }

    @NotNull
    private static String hotspotField(@NotNull SomaticVariant variant) {
        switch (variant.hotspot()) {
            case HOTSPOT:
                return "Yes";
            case NEAR_HOTSPOT:
                return "Near";
            default:
                return Strings.EMPTY;
        }
    }

    @NotNull
    private static String clonalityField(@NotNull EnrichedSomaticVariant variant) {
        switch (variant.clonality()) {
            case CLONAL:
                return "Clonal";
            case SUBCLONAL:
                return "Subclonal";
            case INCONSISTENT:
                return "Inconsistent";
            default:
                return Strings.EMPTY;
        }
    }
}
