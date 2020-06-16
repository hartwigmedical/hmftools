package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.cgi.Cgi;
import com.hartwig.hmftools.vicc.datamodel.civic.Civic;
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
        switch (source()) {
            case ONCOKB: {
                String transcriptId = ((OncoKb) kbSpecificObject()).transcriptID();
                return !transcriptId.isEmpty() ? transcriptId : null;
            }
            case CIVIC: {
                String transcriptId = ((Civic) kbSpecificObject()).coordinates().representativeTranscript();
                if (transcriptId == null) {
                    return null;
                }
                int versionSeparator = transcriptId.indexOf(".");
                return versionSeparator >= 0 ? transcriptId.substring(1, versionSeparator - 1) : transcriptId;
            }
            case CGI: {
                List<String> transcripts = ((Cgi) kbSpecificObject()).transcripts();
                // Even though a vicc entry could have multiple different transcripts, in practice every vicc entry with more than one
                // transcript only has one distinct transcript.
                return !transcripts.isEmpty() ? transcripts.get(0) : null;
            }
            case JAX:
            default:
                return null;
        }
    }
}

