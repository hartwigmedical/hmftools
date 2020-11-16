package com.hartwig.hmftools.serve.sources.docm;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.HotspotFunctions;
import com.hartwig.hmftools.serve.hotspot.ImmutableKnownHotspot;
import com.hartwig.hmftools.serve.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;

import org.jetbrains.annotations.NotNull;

public class DocmExtractor {

    @NotNull
    private final ProteinResolver proteinResolver;

    public DocmExtractor(@NotNull final ProteinResolver proteinResolver) {
        this.proteinResolver = proteinResolver;
    }

    @NotNull
    public List<KnownHotspot> extractFromDocmEntries(@NotNull List<DocmEntry> entries) {
        List<KnownHotspot> knownHotspots = Lists.newArrayList();
        for (DocmEntry entry : entries) {
            List<VariantHotspot> hotspots =
                    proteinResolver.extractHotspotsFromProteinAnnotation(entry.gene(), entry.transcript(), entry.proteinAnnotation());
            for (VariantHotspot hotspot : hotspots) {
                knownHotspots.add(ImmutableKnownHotspot.builder()
                        .from(hotspot)
                        .addSources(Knowledgebase.DOCM)
                        .gene(entry.gene())
                        .transcript(entry.transcript())
                        .proteinAnnotation(entry.proteinAnnotation())
                        .build());
            }
        }

        // Hotspots could appear multiple times in the DoCM. We need to consolidate even though there is only one source.
        return HotspotFunctions.consolidateHotspots(knownHotspots);
    }
}
