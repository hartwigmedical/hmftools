package com.hartwig.hmftools.patientreporter;

import java.io.IOException;

import com.hartwig.hmftools.common.center.Center;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.cosmic.Cosmic;
import com.hartwig.hmftools.common.cosmic.CosmicModel;
import com.hartwig.hmftools.common.ecrf.doid.TumorLocationDoidMapping;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.hmfslicer.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientreporter.data.COSMICGeneFusionModel;
import com.hartwig.hmftools.patientreporter.data.COSMICGeneFusions;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;

import org.jetbrains.annotations.NotNull;

public final class HmfReporterDataLoader {
    private HmfReporterDataLoader() {
    }

    @NotNull
    public static HmfReporterData buildFromFiles(@NotNull final String drupFilterFile, @NotNull final String cosmicFile,
            @NotNull final String centerFile, @NotNull final String signaturePath, @NotNull final String fusionFile)
            throws IOException, HartwigException {
        final GeneModel geneModel = new GeneModel(HmfGenePanelSupplier.asMap());
        final DrupFilter drupFilter = new DrupFilter(drupFilterFile);
        final CosmicModel cosmicModel = Cosmic.buildModelFromCsv(cosmicFile);
        final CenterModel centerModel = Center.readFromCSV(centerFile);
        final COSMICGeneFusionModel fusionModel = COSMICGeneFusions.readFromCSV(fusionFile);
        final TumorLocationDoidMapping doidMapping = TumorLocationDoidMapping.fromResource("/tumor_location_doid_mapping.csv");
        return new ImmutableHmfReporterData(geneModel, cosmicModel, drupFilter, centerModel, signaturePath, fusionModel, doidMapping);
    }
}
