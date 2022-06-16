package com.hartwig.hmftools.orange.algo.protect;

import java.util.List;

import com.hartwig.hmftools.common.protect.ProtectEvidence;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ProtectInterpretedData {

    @NotNull
    public abstract List<ProtectEvidence> reportableEvidences();

    @NotNull
    public abstract List<ProtectEvidence> reportableTrials();

    @NotNull
    public abstract List<ProtectEvidence> unreportedEvidences();

    @NotNull
    public abstract List<ProtectEvidence> unreportedTrials();
}
