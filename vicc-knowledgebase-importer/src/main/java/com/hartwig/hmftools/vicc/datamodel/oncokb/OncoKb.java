package com.hartwig.hmftools.vicc.datamodel.oncokb;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OncoKb implements KbSpecificObject {

    @Nullable
    public abstract OncoKbBiological oncoKbBiological();

    @Nullable
    public abstract OncoKbClinical oncoKbClinical();

    @NotNull
    @Value.Derived
    public String transcriptID() {
        if (oncoKbClinical() != null) {
            return oncoKbClinical().isoform();
        } else if (oncoKbBiological() != null) {
            return oncoKbBiological().isoform();
        } else {
            throw new IllegalStateException("Both biological and clinical record null in OncoKB record!");
        }
    }
}
