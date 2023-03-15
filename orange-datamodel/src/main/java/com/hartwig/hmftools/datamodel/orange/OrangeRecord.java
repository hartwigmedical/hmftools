package com.hartwig.hmftools.datamodel.orange;

import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.time.LocalDate;
import java.util.Set;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeRecord {

    @NotNull
    String sampleId();

    @NotNull
    LocalDate experimentDate();

    @NotNull
    OrangeRefGenomeVersion refGenomeVersion();

    @NotNull
    PurpleRecord purple();

    @NotNull
    LinxRecord linx();

    @NotNull
    LilacRecord lilac();

    @Nullable
    VirusInterpreterData virusInterpreter();

    @Nullable
    ChordRecord chord();

    @Nullable
    CuppaData cuppa();

    @Nullable
    Set<PeachGenotype> peach();

    @NotNull
    OrangePlots plots();

}
