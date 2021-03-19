package com.hartwig.hmftools.serve.refgenome;

import java.util.Map;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.RefGenomeVersion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class RefGenomeResource {

    @NotNull
    public abstract String fastaFile();

    @NotNull
    public abstract Map<String, HmfTranscriptRegion> canonicalTranscriptPerGeneMap();

    @NotNull
    public abstract Map<RefGenomeVersion, String> chainToOtherRefGenomeMap();

}
