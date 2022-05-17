package com.hartwig.hmftools.summon;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.summon.actionability.ActionabilityEntry;
import com.hartwig.hmftools.summon.actionability.ActionabilityKey;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SummonData {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract PurpleData purple();

    @NotNull
    public abstract LinxData linx();

    @NotNull
    public abstract VirusInterpreterData virusInterpreter();

    @NotNull
    public abstract ChordAnalysis chord();

    @NotNull
    public abstract List<ProtectEvidence> protect();

    @NotNull
    public abstract List<ActionabilityEntry> actionabilityEntries();
}