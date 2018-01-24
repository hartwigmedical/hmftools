package com.hartwig.hmftools.patientreporter.report.data;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Comparator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.patientreporter.HmfReporterData;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.FieldBuilder;
import net.sf.dynamicreports.report.datasource.DRDataSource;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.jasperreports.engine.JRDataSource;

public final class GenePanelDataSource {

    public static final FieldBuilder<?> GENE_FIELD = field("gene", String.class);
    public static final FieldBuilder<?> TRANSCRIPT_FIELD = field("transcript", String.class);
    public static final FieldBuilder<?> TYPE_FIELD = field("type", String.class);

    private GenePanelDataSource() {
    }

    @NotNull
    public static JRDataSource fromHmfReporterData(@NotNull final HmfReporterData reporterData) {
        final DRDataSource genePanelDataSource = new DRDataSource(GENE_FIELD.getName(), TRANSCRIPT_FIELD.getName(), TYPE_FIELD.getName());
        final List<HmfGenomeRegion> regions = Lists.newArrayList(reporterData.panelGeneModel().hmfRegions());
        regions.sort(Comparator.comparing(HmfGenomeRegion::gene));

        for (final HmfGenomeRegion region : regions) {
            final String role = reporterData.cosmicGeneModel().getRoleForGene(region.gene());
            genePanelDataSource.add(region.gene(), region.transcript(), role);
        }
        return genePanelDataSource;
    }

    @NotNull
    public static AbstractSimpleExpression<String> transcriptUrl() {
        return new AbstractSimpleExpression<String>() {
            @Override
            public String evaluate(@NotNull final ReportParameters data) {
                return "http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=" + data.getValue(TRANSCRIPT_FIELD.getName());
            }
        };
    }
}
