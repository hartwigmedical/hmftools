package com.hartwig.hmftools.orange.algo;

import java.util.Set;

import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.virusinterpreter.VirusInterpreterData;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OrangeReport {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract Set<DoidNode> configuredTumorLocation();

    @NotNull
    public abstract String cuppaTumorLocation();

    @NotNull
    public abstract PurpleData purpleData();

    @NotNull
    public abstract LinxData linxData();

    @NotNull
    public abstract VirusInterpreterData virusInterpreterData();

    @NotNull
    public abstract ChordAnalysis chordAnalysis();

    @NotNull
    public abstract OrangePlots plots();

}
