package com.hartwig.hmftools.patientreporter;

import java.io.IOException;

import com.hartwig.hmftools.common.center.Center;
import com.hartwig.hmftools.common.center.CenterModel;
import com.hartwig.hmftools.common.cosmic.Cosmic;
import com.hartwig.hmftools.common.cosmic.CosmicModel;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.gene.GeneModel;
import com.hartwig.hmftools.hmfslicer.HmfGeneRegionSupplier;
import com.hartwig.hmftools.patientreporter.filters.DrupFilter;

import org.jetbrains.annotations.NotNull;

public final class HmfReporterDataLoader {
    private HmfReporterDataLoader() {
    }

    @NotNull
    public static HmfReporterData buildFromFiles(@NotNull final String drupFilterFile, @NotNull final String cosmicFile,
            @NotNull final String centerFile, @NotNull final String signaturePath) throws IOException, HartwigException {
        final GeneModel geneModel = new GeneModel(HmfGeneRegionSupplier.asMap());
        final DrupFilter drupFilter = new DrupFilter(drupFilterFile);
        final CosmicModel cosmicModel = Cosmic.buildModelFromCsv(cosmicFile);
        final CenterModel centerModel = Center.readFromCSV(centerFile);
        return new ImmutableHmfReporterData(geneModel, cosmicModel, drupFilter, centerModel, signaturePath);
    }
}
