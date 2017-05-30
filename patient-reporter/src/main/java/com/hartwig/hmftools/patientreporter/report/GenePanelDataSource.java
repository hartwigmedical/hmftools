package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Collection;
import java.util.List;

import com.google.common.collect.Iterables;
import com.hartwig.hmftools.patientreporter.cosmic.CosmicCensus;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

class GenePanelDataSource {

    static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    static final FieldBuilder<?> TRANSCRIPT_FIELD = field("transcript", String.class);
    static final FieldBuilder<?> ROLE_IN_CANCER_FIELD = field("role_in_cancer", String.class);
    static final FieldBuilder<?> GENE2_FIELD = field("gene2", String.class);
    static final FieldBuilder<?> TRANSCRIPT2_FIELD = field("transcript2", String.class);
    static final FieldBuilder<?> ROLE_IN_CANCER2_FIELD = field("role_in_cancer2", String.class);

    private GenePanelDataSource() {
    }

    @NotNull
    static JRDataSource fromCosmic(@NotNull final Collection<String> genes, @NotNull final CosmicCensus cosmicCensus) {
        final DRDataSource genePanelDataSource = new DRDataSource(GENE_FIELD.getName(), TRANSCRIPT_FIELD.getName(),
                ROLE_IN_CANCER_FIELD.getName(), GENE2_FIELD.getName(), TRANSCRIPT2_FIELD.getName(),
                ROLE_IN_CANCER2_FIELD.getName());

        for (final List<String> genePair : Iterables.partition(genes, 2)) {
            final String gene1 = genePair.get(0);
            final String gene2 = genePair.get(1);
            genePanelDataSource.add(gene1, "transcript", "Onco", gene2, "transcript", "Onco");
        }

        return genePanelDataSource;
    }

    @NotNull
    static FieldBuilder<?>[] genePanelFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, TRANSCRIPT_FIELD, ROLE_IN_CANCER_FIELD };
    }

}
