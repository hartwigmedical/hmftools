package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKb;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ViccEntry {

    @NotNull
    public abstract ViccSource source();

    @NotNull
    public abstract List<String> genes();

    @NotNull
    public abstract List<GeneIdentifier> geneIdentifiers();

    @NotNull
    public abstract List<String> featureNames();

    @NotNull
    public abstract List<Feature> features();

    @NotNull
    public abstract Association association();

    @NotNull
    public abstract List<String> tags();

    @NotNull
    public abstract List<String> devTags();

    @NotNull
    public abstract KbSpecificObject kbSpecificObject();

    @Nullable
    @Value.Derived
    public String transcriptId() {
        if (source() == ViccSource.ONCOKB) {
            String transcriptId = ((OncoKb)kbSpecificObject()).transcriptID();
            return !transcriptId.isEmpty() ? transcriptId : null;
        }

        // TODO Implemented for all sources.
        return null;
    }

}

