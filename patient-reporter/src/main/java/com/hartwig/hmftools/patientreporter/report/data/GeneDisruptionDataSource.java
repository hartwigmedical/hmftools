package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.List;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.jasperreports.engine.JRDataSource;

public class GeneDisruptionDataSource {

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> TRANSCRIPT_FIELD = field("transcript", String.class);
    public static final FieldBuilder<?> POSITION_FIELD = field("position", String.class);
    public static final FieldBuilder<?> SV_TYPE_FIELD = field("type", String.class);
    public static final FieldBuilder<?> SV_PARTNER_POSITION_FIELD = field("partner", String.class);
    public static final FieldBuilder<?> SV_ORIENTATION_FIELD = field("orientation", String.class);
    public static final FieldBuilder<?> SV_GENE_CONTEXT = field("gene context", String.class);
    public static final FieldBuilder<?> SV_VAF = field("vaf", String.class);

    private GeneDisruptionDataSource() {
    }

    @NotNull
    public static JRDataSource fromGeneDisruptions(@NotNull List<GeneDisruptionData> disruptions) {

        final DRDataSource dataSource =
                new DRDataSource(GENE_FIELD.getName(), TRANSCRIPT_FIELD.getName(), POSITION_FIELD.getName(), SV_TYPE_FIELD.getName(),
                        SV_PARTNER_POSITION_FIELD.getName(), SV_ORIENTATION_FIELD.getName(), SV_GENE_CONTEXT.getName(), SV_VAF.getName());

        disruptions.forEach(
                g -> dataSource.add(g.geneName(), g.transcript(), g.location(), g.type(), g.partner(), g.orientation(), g.geneContext(),
                        g.vaf()));

        return dataSource;
    }

    @NotNull
    public static FieldBuilder<?>[] geneDisruptionFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, TRANSCRIPT_FIELD, POSITION_FIELD, SV_TYPE_FIELD, SV_PARTNER_POSITION_FIELD,
                SV_ORIENTATION_FIELD, SV_GENE_CONTEXT, SV_VAF };
    }
}
