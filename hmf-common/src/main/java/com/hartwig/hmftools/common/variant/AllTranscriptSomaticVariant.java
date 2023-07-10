package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

public interface AllTranscriptSomaticVariant extends SomaticVariant {
    List<VariantTranscriptImpact> allTranscripts();
}
