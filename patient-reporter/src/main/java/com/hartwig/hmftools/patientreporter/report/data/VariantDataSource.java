package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public class VariantDataSource {

    private static final String COSMIC_IDENTIFIER = "COSM";

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> POSITION_FIELD = field("position", String.class);
    public static final FieldBuilder<?> VARIANT_FIELD = field("variant", String.class);
    public static final FieldBuilder<?> HGVS_CODING_FIELD = field("hgvs_coding", String.class);
    public static final FieldBuilder<?> HGVS_PROTEIN_FIELD = field("hgvs_protein", String.class);
    public static final FieldBuilder<?> CONSEQUENCE_FIELD = field("consequence", String.class);
    public static final FieldBuilder<?> COSMIC_FIELD = field("cosmic", String.class);
    public static final FieldBuilder<?> COSMIC_NR_FIELD = field("cosmic_nr", String.class);
    public static final FieldBuilder<?> DEPTH_VAF_FIELD = field("depth_vaf", String.class);
    public static final FieldBuilder<?> PLOIDY_TAF_FIELD = field("ploidy_taf", String.class);

    private VariantDataSource() {
    }

    @NotNull
    public static JRDataSource fromVariants(@NotNull final List<VariantReport> variantReports,
            @NotNull final HmfReporterData reporterData) {
        final DRDataSource variantDataSource =
                new DRDataSource(GENE_FIELD.getName(), POSITION_FIELD.getName(), VARIANT_FIELD.getName(), DEPTH_VAF_FIELD.getName(),
                        COSMIC_FIELD.getName(), COSMIC_NR_FIELD.getName(), HGVS_CODING_FIELD.getName(), HGVS_PROTEIN_FIELD.getName(),
                        CONSEQUENCE_FIELD.getName(), PLOIDY_TAF_FIELD.getName());

        for (final VariantReport variantReport : variantReports) {
            final String displayGene = reporterData.drupFilter().test(variantReport) ? variantReport.gene() + " *" : variantReport.gene();
            variantDataSource.add(displayGene, variantReport.variant().chromosomePosition(), variantReport.variantField(),
                    variantReport.depthVafField(), variantReport.cosmicID(), stripCosmicIdentifier(variantReport.cosmicID()),
                    variantReport.hgvsCoding(), variantReport.hgvsProtein(), variantReport.consequence(), variantReport.ploidyTafField());
        }

        return variantDataSource;
    }

    @NotNull
    private static String stripCosmicIdentifier(@NotNull final String cosmicID) {
        final int identifierPos = cosmicID.indexOf(COSMIC_IDENTIFIER);
        if (identifierPos >= 0) {
            return cosmicID.substring(identifierPos + COSMIC_IDENTIFIER.length());
        } else {
            return cosmicID;
        }
    }

    @NotNull
    public static FieldBuilder<?>[] variantFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, POSITION_FIELD, VARIANT_FIELD, HGVS_CODING_FIELD, HGVS_PROTEIN_FIELD, CONSEQUENCE_FIELD,
                COSMIC_FIELD, COSMIC_NR_FIELD, DEPTH_VAF_FIELD, PLOIDY_TAF_FIELD };
    }
}
