package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.driver.Driver;
import com.hartwig.hmftools.datamodel.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface Virus extends Driver {

    @NotNull
    VirusInterpreterEntry interpreterEntry();

    @NotNull
    default String name()
    {
        return interpreterEntry().name();
    }

    @NotNull
    default VirusBreakendQCStatus qcStatus()
    {
        return interpreterEntry().qcStatus();
    }

    default int integrations()
    {
        return interpreterEntry().integrations();
    }

    @Nullable
    default VirusInterpretation interpretation()
    {
        return interpreterEntry().interpretation();
    }

    default double percentageCovered()
    {
        return interpreterEntry().percentageCovered();
    }

    default double meanCoverage()
    {
        return interpreterEntry().meanCoverage();
    }

    @Nullable
    default Double expectedClonalCoverage()
    {
        return interpreterEntry().expectedClonalCoverage();
    }
}
