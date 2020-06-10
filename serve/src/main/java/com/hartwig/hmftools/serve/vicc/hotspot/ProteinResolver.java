package com.hartwig.hmftools.serve.vicc.hotspot;

import java.util.List;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface ProteinResolver {

    @NotNull
    List<VariantHotspot> extractHotspotsFromProteinAnnotation(@NotNull String gene, @Nullable String specificTranscript,
            @NotNull String proteinAnnotation);
}
