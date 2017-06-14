package com.hartwig.hmftools.patientreporter.report;

import static com.google.common.base.Strings.nullToEmpty;

import static net.sf.dynamicreports.report.builder.DynamicReports.field;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.patientreporter.HmfReporterData;

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
    static JRDataSource fromHmfReporterData(@NotNull final HmfReporterData reporterData) {
        final DRDataSource genePanelDataSource = new DRDataSource(GENE_FIELD.getName(), TRANSCRIPT_FIELD.getName(),
                TYPE_FIELD.getName(), GENE2_FIELD.getName(), TRANSCRIPT2_FIELD.getName(), TYPE2_FIELD.getName());
        final List<HmfGenomeRegion> regions = reporterData.slicer().hmfRegions().stream().collect(Collectors.toList());
        regions.sort(Comparator.comparing(HmfGenomeRegion::gene));

        for (int i = 0; i < regions.size() / 2; i++) {
            final HmfGenomeRegion region1 = regions.get(i);
            final String role1 = reporterData.cosmicModel().getRoleForGene(region1.gene());
            final int secondIndex = regions.size() / 2 + i;
            final HmfGenomeRegion region2 = secondIndex < regions.size() ? regions.get(secondIndex) : null;

            final String geneName2 = region2 != null ? region2.gene() : Strings.EMPTY;
            final String transcript2 = region2 != null ? region2.transcript() : Strings.EMPTY;
            final String role2 =
                    region2 != null ? reporterData.cosmicModel().getRoleForGene(region2.gene()) : Strings.EMPTY;

            genePanelDataSource.add(region1.gene(), region1.transcript(), role1, geneName2, transcript2,
                    role2);
        }

        return genePanelDataSource;
    }

    @NotNull
    static FieldBuilder<?>[] genePanelFields() {
        return new FieldBuilder<?>[] { GENE_FIELD, TRANSCRIPT_FIELD, TYPE_FIELD, GENE2_FIELD, TRANSCRIPT2_FIELD,
                TYPE2_FIELD };
    }
}
