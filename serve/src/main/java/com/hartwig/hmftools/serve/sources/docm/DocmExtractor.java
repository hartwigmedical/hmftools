package com.hartwig.hmftools.serve.sources.docm;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;

import org.jetbrains.annotations.NotNull;

public class DocmExtractor {

    @NotNull
    private final ProteinResolver proteinResolver;

    public DocmExtractor(@NotNull final ProteinResolver proteinResolver) {
        this.proteinResolver = proteinResolver;
    }

    @NotNull
    public Map<DocmEntry, List<VariantHotspot>> extractFromDocmEntries(@NotNull List<DocmEntry> entries) {
        Map<DocmEntry, List<VariantHotspot>> hotspotsPerEntry = Maps.newHashMap();
        for (DocmEntry entry : entries) {
            hotspotsPerEntry.put(entry,
                    proteinResolver.extractHotspotsFromProteinAnnotation(entry.gene(), entry.transcript(), entry.proteinAnnotation()));
        }
        return hotspotsPerEntry;
    }
}
