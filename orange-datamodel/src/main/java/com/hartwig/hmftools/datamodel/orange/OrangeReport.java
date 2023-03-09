package com.hartwig.hmftools.datamodel.orange;

import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.time.LocalDate;
import java.util.List;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OrangeReport {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract LocalDate experimentDate();

    @NotNull
    public abstract OrangeRefGenomeVersion refGenomeVersion();

    @NotNull
    public abstract PurpleRecord purple();

    @NotNull
    public abstract LinxRecord linx();

    @NotNull
    public abstract LilacRecord lilac();

    @Nullable
    public abstract VirusInterpreterData virusInterpreter();

    @Nullable
    public abstract ChordRecord chord();

    @Nullable
    public abstract CuppaData cuppa();

    @Nullable
    public abstract List<PeachGenotype> peach();

    @NotNull
    public abstract OrangePlots plots();

}
