package com.hartwig.hmftools.knowledgebasegenerator.hotspot;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.knowledgebasegenerator.RefGenomeVersion;
import com.hartwig.hmftools.knowledgebasegenerator.transvar.Transvar;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class HotspotExtractor {

    @NotNull
    private final Transvar transvar;

    public static HotspotExtractor fromRefGenome(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile) {
        return new HotspotExtractor(new Transvar(refGenomeVersion, refGenomeFastaFile));
    }

    private HotspotExtractor(@NotNull Transvar transvar) {
        this.transvar = transvar;
    }

    @NotNull
    public List<VariantHotspot> extractHotspots(@NotNull ViccEntry viccEntry) {
        // TODO Implement

        return Lists.newArrayList();
    }
}
