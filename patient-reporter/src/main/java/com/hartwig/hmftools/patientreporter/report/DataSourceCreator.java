package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

class DataSourceCreator {

    static final String GENE_FIELD = "gene";
    static final String POSITION_FIELD = "position";
    static final String VARIANT_FIELD = "variant";
    static final String TRANSCRIPT_FIELD = "transcript";
    static final String HGVS_CODING_FIELD = "hgvs_coding";
    static final String HGVS_PROTEIN_FIELD = "hgvs_protein";
    static final String EFFECT_FIELD = "effect";
    static final String COSMIC_FIELD = "cosmic";
    static final String READ_COUNT_FIELD = "frequency";

    private DataSourceCreator() {
    }

    @NotNull
    static JRDataSource fromVariants(@NotNull final List<VariantReport> variants) {
        final DRDataSource variantDataSource = new DRDataSource(GENE_FIELD, POSITION_FIELD, VARIANT_FIELD,
                TRANSCRIPT_FIELD, HGVS_CODING_FIELD, HGVS_PROTEIN_FIELD, EFFECT_FIELD, COSMIC_FIELD, READ_COUNT_FIELD);

        for (final VariantReport variant : variants) {
            variantDataSource.add(variant.gene(), variant.position(), variant.ref() + " > " + variant.alt(),
                    variant.transcript(), variant.hgvsCoding(), variant.hgvsProtein(), variant.consequence(),
                    variant.cosmicID(), variant.alleleReadCount() + " / " + variant.totalReadCount());
        }

        return variantDataSource;
    }

    @NotNull
    static FieldBuilder<?>[] variantFields() {
        final FieldBuilder<?>[] fieldList = new FieldBuilder<?>[9];
        fieldList[0] = field(GENE_FIELD, String.class);
        fieldList[1] = field(POSITION_FIELD, String.class);
        fieldList[2] = field(VARIANT_FIELD, String.class);
        fieldList[3] = field(TRANSCRIPT_FIELD, String.class);
        fieldList[4] = field(HGVS_CODING_FIELD, String.class);
        fieldList[5] = field(HGVS_PROTEIN_FIELD, String.class);
        fieldList[6] = field(EFFECT_FIELD, String.class);
        fieldList[7] = field(COSMIC_FIELD, String.class);
        fieldList[8] = field(READ_COUNT_FIELD, String.class);
        return fieldList;
    }
}
