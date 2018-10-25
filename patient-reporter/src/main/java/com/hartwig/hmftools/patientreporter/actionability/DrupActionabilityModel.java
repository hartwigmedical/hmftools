package com.hartwig.hmftools.patientreporter.actionability;

import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DrupActionabilityModel {

    @NotNull
    public abstract Set<String> actionableGenes();
    @NotNull
    public abstract Map<String, DriverCategory> geneDriverCategoryMap();

}
