package com.hartwig.hmftools.patientdb.data;

import java.util.List;

import com.hartwig.hmftools.common.doid.DoidNode;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class CuratedPrimaryTumor {

    @Nullable
    public abstract String searchTerm();

    @Nullable
    public abstract String location();

    @Nullable
    public abstract String subLocation();

    @Nullable
    public abstract String type();

    @Nullable
    public abstract String subType();

    @Nullable
    public abstract String extraDetails();

    @Nullable
    public abstract List<DoidNode> doidNodes();

    @Nullable
    public abstract List<String> snomedIds();

}
