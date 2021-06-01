package com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class FormStatus {

    @NotNull
    public abstract FormStatusState state();

    public abstract boolean locked();

    @NotNull
    @Value.Derived
    public String stateString() {
        return state().stateString();
    }

    @Nullable
    @Value.Derived
    public String lockedString() {
        return state().equals(FormStatusState.UNDEFINED) ? null : Boolean.toString(locked());
    }

    @NotNull
    public static FormStatus undefined() {
        return ImmutableFormStatus.builder().state(FormStatusState.UNDEFINED).locked(false).build();
    }

    @NotNull
    public static FormStatus merge(@NotNull FormStatus... states) {
        assert states.length >= 1;

        FormStatusState bestState = states[0].state();
        boolean mergedLocked = states[0].locked();
        for (int i = 1; i < states.length; i++) {
            bestState = FormStatusState.best(bestState, states[i].state());
            mergedLocked = mergedLocked && states[i].locked();
        }
        return ImmutableFormStatus.builder().state(bestState).locked(mergedLocked).build();
    }
}
