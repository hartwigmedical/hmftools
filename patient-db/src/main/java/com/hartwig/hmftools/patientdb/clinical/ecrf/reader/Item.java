package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class Item implements OIDObject {

    @Override
    @Value.Parameter
    @NotNull
    public abstract String oid();

    @Override
    @Value.Parameter
    @NotNull
    public abstract String name();

    @Value.Parameter
    @Nullable
    public abstract String codeListOID();
}
