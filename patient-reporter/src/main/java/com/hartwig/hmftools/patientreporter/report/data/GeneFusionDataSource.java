package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public final class GeneFusionDataSource {

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> TRANSCRIPT_FIELD = field("transcript", String.class);
    public static final FieldBuilder<?> POSITION_FIELD = field("position", String.class);
    public static final FieldBuilder<?> SV_TYPE_FIELD = field("type", String.class);
    public static final FieldBuilder<?> SV_PARTNER_GENE_FIELD = field("partner_gene", String.class);
    public static final FieldBuilder<?> SV_PARTNER_TRANSCRIPT_FIELD = field("partner_transcript", String.class);
    public static final FieldBuilder<?> SV_PARTNER_POSITION_FIELD = field("partner", String.class);
    public static final FieldBuilder<?> SV_PARTNER_CONTEXT_FIELD = field("partner_context", String.class);
    public static final FieldBuilder<?> SV_GENE_CONTEXT = field("gene context", String.class);
    public static final FieldBuilder<?> SV_VAF = field("vaf", String.class);

    private GeneFusionDataSource() {
    }

    @NotNull
    public static JRDataSource fromGeneFusions(@NotNull List<GeneFusionData> fusions) {
        final DRDataSource dataSource =
                new DRDataSource(GENE_FIELD.getName(), TRANSCRIPT_FIELD.getName(), POSITION_FIELD.getName(), SV_GENE_CONTEXT.getName(),
                        SV_PARTNER_GENE_FIELD.getName(), SV_PARTNER_TRANSCRIPT_FIELD.getName(), SV_PARTNER_POSITION_FIELD.getName(),
                        SV_PARTNER_CONTEXT_FIELD.getName(), SV_TYPE_FIELD.getName(), SV_VAF.getName());

        fusions.forEach(
                g -> dataSource.add(g.geneStart(), g.transcriptStart(), g.start(), g.geneContextStart(), g.geneEnd(), g.transcriptEnd(),
                        g.end(), g.geneContextEnd(), g.type(), g.vaf()));

        return dataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] geneFusionFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, TRANSCRIPT_FIELD, POSITION_FIELD, SV_GENE_CONTEXT, SV_PARTNER_GENE_FIELD,
                SV_PARTNER_TRANSCRIPT_FIELD, SV_PARTNER_POSITION_FIELD, SV_PARTNER_CONTEXT_FIELD, SV_TYPE_FIELD, SV_VAF };
    }
}
