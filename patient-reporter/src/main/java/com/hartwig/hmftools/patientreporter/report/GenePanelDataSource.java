package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Collection;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.patientreporter.genePanel.GenePanelModel;
import com.hartwig.hmftools.patientreporter.slicing.HMFSlicingAnnotation;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

class GenePanelDataSource {

    static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    static final FieldBuilder<?> TRANSCRIPT_FIELD = field("transcript", String.class);
    static final FieldBuilder<?> GENE_TYPE_FIELD = field("gene_type", String.class);
    static final FieldBuilder<?> GENE2_FIELD = field("gene2", String.class);
    static final FieldBuilder<?> TRANSCRIPT2_FIELD = field("transcript2", String.class);
    static final FieldBuilder<?> GENE2_TYPE_FIELD = field("gene2_type", String.class);

    private GenePanelDataSource() {
    }

    @NotNull
    static JRDataSource fromCosmic(@NotNull final Collection<HMFSlicingAnnotation> genes,
            @NotNull final GenePanelModel genePanelModel) {
        final DRDataSource genePanelDataSource = new DRDataSource(GENE_FIELD.getName(), TRANSCRIPT_FIELD.getName(),
                GENE_TYPE_FIELD.getName(), GENE2_FIELD.getName(), TRANSCRIPT2_FIELD.getName(),
                GENE2_TYPE_FIELD.getName());

        for (final List<HMFSlicingAnnotation> annotationPair : Iterables.partition(genes, 2)) {
            final HMFSlicingAnnotation annotation1 = annotationPair.get(0);
            final HMFSlicingAnnotation annotation2 = annotationPair.get(1);
            genePanelDataSource.add(annotation1.gene(), annotation1.transcript(),
                    genePanelModel.type(annotation1.gene()), annotation2.gene(), annotation2.transcript(),
                    genePanelModel.type(annotation2.gene()));
        }

        return genePanelDataSource;
    }

    @NotNull
    static FieldBuilder<?>[] genePanelFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, TRANSCRIPT_FIELD, GENE_TYPE_FIELD };
    }

}
