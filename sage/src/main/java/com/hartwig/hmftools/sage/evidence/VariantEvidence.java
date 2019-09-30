package com.hartwig.hmftools.sage.evidence;

import java.util.List;

import com.hartwig.hmftools.common.hotspot.VariantHotspot;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface VariantEvidence extends VariantHotspot {

    // INFO FIELDS
    @NotNull
    String readContext();

    // FORMAT FIELDS
    @NotNull
    SampleEvidence normalEvidence();

    @NotNull
    List<SampleEvidence> tumorEvidence();

}
