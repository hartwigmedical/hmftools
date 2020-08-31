package com.hartwig.hmftools.serve.docm;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinToHotspotConverter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class DocmExtractor {

    private static final Logger LOGGER = LogManager.getLogger(DocmExtractor.class);

    @NotNull
    private final ProteinToHotspotConverter proteinToHotspotConverter;

    public DocmExtractor(@NotNull final ProteinToHotspotConverter proteinToHotspotConverter) {
        this.proteinToHotspotConverter = proteinToHotspotConverter;
    }

    @NotNull
    public Map<DocmEntry, List<VariantHotspot>> extractFromDocmEntries(@NotNull List<DocmEntry> entries) {
        Map<DocmEntry, List<VariantHotspot>> hotspotsPerEntry = Maps.newHashMap();
        for (DocmEntry entry : entries) {
            if (ProteinToHotspotConverter.isResolvableProteinAnnotation(entry.proteinAnnotation())) {
                hotspotsPerEntry.put(entry, proteinToHotspotConverter.resolveProteinAnnotation(entry.gene(), entry.transcript(), entry.proteinAnnotation()));
            } else {
                LOGGER.warn("Cannot resolve DoCM protein annotation: '{}:p.{}'", entry.gene(), entry.proteinAnnotation());
            }
        }
        return hotspotsPerEntry;
    }
}
