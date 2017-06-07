package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

class GenePanelDataSource {

    static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    static final FieldBuilder<?> TRANSCRIPT_FIELD = field("transcript", String.class);
    static final FieldBuilder<?> TYPE_FIELD = field("type", String.class);
    static final FieldBuilder<?> GENE2_FIELD = field("gene2", String.class);
    static final FieldBuilder<?> TRANSCRIPT2_FIELD = field("transcript2", String.class);
    static final FieldBuilder<?> TYPE2_FIELD = field("type2", String.class);

    private GenePanelDataSource() {
    }

    @NotNull
    static JRDataSource fromHMFSlicingAnnotations(@NotNull final List<HMFSlicingAnnotation> genes) {
        final DRDataSource genePanelDataSource = new DRDataSource(GENE_FIELD.getName(), TRANSCRIPT_FIELD.getName(),
                TYPE_FIELD.getName(), GENE2_FIELD.getName(), TRANSCRIPT2_FIELD.getName(), TYPE2_FIELD.getName());

        for (int i = 0; i < genes.size() / 2; i++) {
            final HMFSlicingAnnotation gene1 = genes.get(i);
            final int secondIndex = genes.size() / 2 + i;
            final HMFSlicingAnnotation gene2 = secondIndex < genes.size() ? genes.get(secondIndex) : null;

            final String geneName2 = gene2 != null ? gene2.gene() : Strings.EMPTY;
            final String transcript2 = gene2 != null ? gene2.transcript() : Strings.EMPTY;

            genePanelDataSource.add(gene1.gene(), gene1.transcript(), Strings.EMPTY, geneName2, transcript2,
                    Strings.EMPTY);
        }

        return genePanelDataSource;
    }

    @NotNull
    static FieldBuilder<?>[] genePanelFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, TRANSCRIPT_FIELD, TYPE_FIELD, GENE2_FIELD, TRANSCRIPT2_FIELD,
                TYPE2_FIELD };
    }
}
