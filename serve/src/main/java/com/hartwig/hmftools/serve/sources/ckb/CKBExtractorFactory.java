package com.hartwig.hmftools.serve.sources.ckb;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.curation.DoidLookup;
import com.hartwig.hmftools.serve.extraction.EventExtractorFactory;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResource;
import com.hartwig.hmftools.serve.sources.vicc.curation.DrugCurator;
import com.hartwig.hmftools.serve.sources.vicc.curation.EvidenceLevelCurator;

import org.jetbrains.annotations.NotNull;

public class CKBExtractorFactory {

    private CKBExtractorFactory() {

    }

    @NotNull
    public static CKBExtractor buildCkbExtractor(@NotNull EventClassifierConfig config, @NotNull RefGenomeResource refGenomeResource,
            @NotNull List<DriverGene> driverGenes, @NotNull KnownFusionCache knownFusionCache, @NotNull DoidLookup missingDoidLookup) {
        return new CKBExtractor(EventExtractorFactory.create(config, refGenomeResource, driverGenes, knownFusionCache),
                new ActionableEvidenceFactory(missingDoidLookup, new DrugCurator(), new EvidenceLevelCurator()));
    }
}
