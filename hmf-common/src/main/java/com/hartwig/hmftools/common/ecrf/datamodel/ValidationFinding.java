package com.hartwig.hmftools.common.ecrf.datamodel;

import org.apache.logging.log4j.message.Message;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ValidationFinding implements Message {
    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String ecrfItem();

    @NotNull
    public abstract String message();

    @NotNull
    public abstract String formStatus();

    @NotNull
    public abstract String formLocked();

    @Override
    public String getFormat() {
        return message();
    }

    @Override
    public String getFormattedMessage() {
        return message();
    }

    @Override
    public Object[] getParameters() {
        return null;
    }

    @Override
    public Throwable getThrowable() {
        return null;
    }
}
