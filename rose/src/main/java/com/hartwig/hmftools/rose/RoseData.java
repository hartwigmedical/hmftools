package com.hartwig.hmftools.rose;

import java.util.List;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.clinical.PatientPrimaryTumor;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOrigin;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.purple.PurpleData;
import com.hartwig.hmftools.common.virus.VirusInterpreterData;
import com.hartwig.hmftools.rose.actionability.ActionabilityEntry;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class RoseData {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract PurpleData purple();

    @NotNull
    public abstract LinxData linx();

    @NotNull
    public abstract VirusInterpreterData virusInterpreter();

    @NotNull
    public abstract ChordData chord();

    @NotNull
    public abstract MolecularTissueOrigin molecularTissueOrigin();

    @NotNull
    public abstract List<ActionabilityEntry> actionabilityEntries();

    @NotNull
    public abstract List<DriverGene> driverGenes();
    @NotNull
    public abstract List<PatientPrimaryTumor> patientPrimaryTumors();
}