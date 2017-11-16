package com.hartwig.hmftools.patientreporter;

import java.io.IOException;

import com.hartwig.hmftools.common.cosmic.Cosmic;
import com.hartwig.hmftools.common.cosmic.CosmicModel;
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
            @NotNull final String fusionFile) throws IOException, HartwigException {
        final GeneModel geneModel = new GeneModel(HmfGenePanelSupplier.hmfGeneMap());
        final DrupFilter drupFilter = new DrupFilter(drupFilterFile);
        final CosmicModel cosmicModel = Cosmic.buildModelFromCsv(cosmicFile);
        final COSMICGeneFusionModel fusionModel = COSMICGeneFusions.readFromCSV(fusionFile);
        return ImmutableHmfReporterData.of(geneModel, cosmicModel, drupFilter, fusionModel);
    }
}
