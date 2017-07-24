package com.hartwig.hmftools.common.ecrf.formstatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class FormStatusData {

    @NotNull
    public abstract String dataStatus();

    @NotNull
    public abstract String locked();

    @NotNull
    public String dataStatusString() {
        switch (dataStatus()) {
            case "0":
                return "saved";
            case "1":
                return "submitted";
            case "2":
                return "submitted with discrepancies";
            case "3":
                return "submitted with missing";
            case "4":
                return "verified";
            default:
                return "unknown";
        }
    }
}
