package com.hartwig.hmftools.knowledgebasegenerator.hotspot;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.knowledgebasegenerator.RefGenomeVersion;
import com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase.Source;
import com.hartwig.hmftools.knowledgebasegenerator.transvar.Transvar;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class HotspotExtractor {

    private static final Logger LOGGER = LogManager.getLogger(HotspotExtractor.class);

    @NotNull
    private final Transvar transvar;

    public static HotspotExtractor fromRefGenome(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile) {
        return new HotspotExtractor(new Transvar(refGenomeVersion, refGenomeFastaFile));
    }

    private HotspotExtractor(@NotNull Transvar transvar) {
        this.transvar = transvar;
    }

    @NotNull
    public List<VariantHotspot> extractHotspots(@NotNull ViccEntry viccEntry) throws IOException, InterruptedException {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        if (Source.sourceFromKnowledgebase(viccEntry.source()) == Source.ONCOKB) {
            for (Feature feature : viccEntry.features()) {
                if (feature.biomarkerType().equals("missense_variant")) {
                    LOGGER.info("Converting feature on {} with name '{}'", feature.geneSymbol(), feature.name());
                    hotspots.addAll(transvar.extractHotspotsFromProteinAnnotation(feature.geneSymbol(), feature.name()));
                } else {
                    LOGGER.info("Skipping feature interpretation because of biomarket type '{}'", feature.biomarkerType());
                }
            }
        }

        return hotspots;
    }
}
