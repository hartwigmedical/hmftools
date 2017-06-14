package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

class PatientDataSource {

    private static final String COSMIC_IDENTIFIER = "COSM";

    static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);

    static final FieldBuilder<?> POSITION_FIELD = field("position", String.class);
    static final FieldBuilder<?> VARIANT_FIELD = field("variant", String.class);
    static final FieldBuilder<?> HGVS_CODING_FIELD = field("hgvs_coding", String.class);
    static final FieldBuilder<?> HGVS_PROTEIN_FIELD = field("hgvs_protein", String.class);
    static final FieldBuilder<?> CONSEQUENCE_FIELD = field("consequence", String.class);
    static final FieldBuilder<?> COSMIC_FIELD = field("cosmic", String.class);
    static final FieldBuilder<?> COSMIC_NR_FIELD = field("cosmic_nr", String.class);
    static final FieldBuilder<?> DEPTH_VAF_FIELD = field("depth_vaf", String.class);
    static final FieldBuilder<?> PLOIDY_TAF_FIELD = field("ploidy_taf", String.class);

    static final FieldBuilder<?> CHROMOSOME_FIELD = field("chromosome", String.class);
    static final FieldBuilder<?> BAND_FIELD = field("band", String.class);
    static final FieldBuilder<?> COPY_NUMBER_TYPE_FIELD = field("copynumber_type", String.class);
    static final FieldBuilder<?> COPY_NUMBER_FIELD = field("copynumber", String.class);

    private PatientDataSource() {
    }

    @NotNull
    static JRDataSource fromVariants(@NotNull final List<VariantReport> variants,
            @NotNull final DrupFilter drupFilter) {
        final DRDataSource variantDataSource = new DRDataSource(GENE_FIELD.getName(), POSITION_FIELD.getName(),
                VARIANT_FIELD.getName(), DEPTH_VAF_FIELD.getName(), COSMIC_FIELD.getName(), COSMIC_NR_FIELD.getName(),
                HGVS_CODING_FIELD.getName(), HGVS_PROTEIN_FIELD.getName(), CONSEQUENCE_FIELD.getName(),
                PLOIDY_TAF_FIELD.getName());

        for (final VariantReport variant : variants) {
            final String displayGene = drupFilter.test(variant) ? variant.gene() + " *" : variant.gene();
            variantDataSource.add(displayGene, variant.chromosomePosition(), variant.variantField(),
                    variant.depthVafField(), variant.cosmicID(), stripCosmicIdentifier(variant.cosmicID()),
                    variant.hgvsCoding(), variant.hgvsProtein(), variant.consequence(), variant.ploidyTafField());
        }

        return variantDataSource;
    }

    @NotNull
    static JRDataSource fromCopyNumbers(@NotNull final List<CopyNumberReport> copyNumbers) {
        final DRDataSource copyNumberDatasource = new DRDataSource(CHROMOSOME_FIELD.getName(), BAND_FIELD.getName(),
                GENE_FIELD.getName(), COPY_NUMBER_TYPE_FIELD.getName(), COPY_NUMBER_FIELD.getName());

        for (final CopyNumberReport copyNumber : copyNumbers) {
            copyNumberDatasource.add(copyNumber.chromosome(), copyNumber.chromosomeBand(), copyNumber.gene(),
                    copyNumber.resolveType(), Integer.toString(copyNumber.copyNumber()));
        }
        return copyNumberDatasource;
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
    static FieldBuilder<?>[] variantFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, POSITION_FIELD, VARIANT_FIELD, HGVS_CODING_FIELD,
                HGVS_PROTEIN_FIELD, CONSEQUENCE_FIELD, COSMIC_FIELD, COSMIC_NR_FIELD, DEPTH_VAF_FIELD,
                PLOIDY_TAF_FIELD };
    }

    @NotNull
    static FieldBuilder<?>[] copyNumberFields() {
        return new FieldBuilder<?>[] { CHROMOSOME_FIELD, BAND_FIELD, GENE_FIELD, COPY_NUMBER_TYPE_FIELD,
                COPY_NUMBER_FIELD };
    }
}
