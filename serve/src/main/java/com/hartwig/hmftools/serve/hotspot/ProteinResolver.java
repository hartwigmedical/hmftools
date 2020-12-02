package com.hartwig.hmftools.serve.hotspot;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface ProteinResolver {

    @NotNull
    List<VariantHotspot> resolve(@NotNull String gene, @Nullable String specificTranscript, @NotNull String proteinAnnotation);

    @NotNull
    Set<String> unresolvedProteinAnnotations();
}
