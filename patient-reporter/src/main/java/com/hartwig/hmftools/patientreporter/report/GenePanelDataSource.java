package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Collection;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

class GenePanelDataSource {

    static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    static final FieldBuilder<?> TRANSCRIPT_FIELD = field("transcript", String.class);
    static final FieldBuilder<?> GENE2_FIELD = field("gene2", String.class);
    static final FieldBuilder<?> TRANSCRIPT2_FIELD = field("transcript2", String.class);
    static final FieldBuilder<?> GENE3_FIELD = field("gene3", String.class);
    static final FieldBuilder<?> TRANSCRIPT3_FIELD = field("transcript3", String.class);

    private GenePanelDataSource() {
    }

    @NotNull
    static JRDataSource fromHMFSlicingAnnotations(@NotNull final Collection<HMFSlicingAnnotation> genes) {
        final DRDataSource genePanelDataSource = new DRDataSource(GENE_FIELD.getName(), TRANSCRIPT_FIELD.getName(),
                GENE2_FIELD.getName(), TRANSCRIPT2_FIELD.getName(), GENE3_FIELD.getName(),
                TRANSCRIPT3_FIELD.getName());

        for (final List<HMFSlicingAnnotation> annotationList : Iterables.paddedPartition(genes, 3)) {
            final HMFSlicingAnnotation annotation1 = annotationList.get(0);
            final HMFSlicingAnnotation annotation2 = annotationList.get(1);
            final HMFSlicingAnnotation annotation3 = annotationList.get(2);

            final String gene2 = annotation2 != null ? annotation2.gene() : Strings.EMPTY;
            final String transcript2 = annotation2 != null ? annotation2.transcript() : Strings.EMPTY;

            final String gene3 = annotation3 != null ? annotation3.gene() : Strings.EMPTY;
            final String transcript3 = annotation3 != null ? annotation3.transcript() : Strings.EMPTY;

            genePanelDataSource.add(annotation1.gene(), annotation1.transcript(), gene2, transcript2, gene3,
                    transcript3);
        }

        return genePanelDataSource;
    }

    @NotNull
    static FieldBuilder<?>[] genePanelFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, TRANSCRIPT_FIELD, GENE2_FIELD, TRANSCRIPT2_FIELD, GENE3_FIELD,
                TRANSCRIPT3_FIELD };
    }
}
