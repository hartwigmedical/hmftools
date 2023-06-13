package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
import org.jetbrains.annotations.NotNull;

import java.util.List;

public interface MultipleTranscriptSomaticVariant extends SomaticVariant {

    @NotNull
    List<VariantTranscriptImpact> transcripts();
}
