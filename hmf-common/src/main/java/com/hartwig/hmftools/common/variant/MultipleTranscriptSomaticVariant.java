package com.hartwig.hmftools.common.variant;

import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
import org.jetbrains.annotations.Nullable;

import java.util.List;

public interface MultipleTranscriptSomaticVariant extends SomaticVariant {

    @Nullable
    List<VariantTranscriptImpact> transcripts();
}
