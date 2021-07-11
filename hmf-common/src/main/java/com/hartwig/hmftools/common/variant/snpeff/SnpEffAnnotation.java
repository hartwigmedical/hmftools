package com.hartwig.hmftools.common.variant.snpeff;

import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.variant.VariantConsequence;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class SnpEffAnnotation
{
    private static final String FEATURE_TYPE_TRANSCRIPT = "transcript";

    @NotNull
    public abstract String gene();

    @NotNull
    public abstract String allele();

    @NotNull
    public abstract String effects();

    @NotNull
    public abstract List<VariantConsequence> consequences();

    @NotNull
    public String consequenceString()
    {
        final StringJoiner consequenceString = new StringJoiner("; ");
        for(final VariantConsequence consequence : consequences())
        {
            if(!consequence.readableSequenceOntologyTerm().isEmpty())
            {
                consequenceString.add(consequence.readableSequenceOntologyTerm());
            }
        }
        return consequenceString.toString();
    }

    @NotNull
    public abstract String severity();

    @NotNull
    public abstract String geneID();

    @NotNull
    public abstract String featureType();

    @NotNull
    public abstract String featureID();

    @NotNull
    abstract String transcriptBioType();

    @NotNull
    public abstract String rank();

    @NotNull
    public abstract String hgvsCoding();

    @NotNull
    public abstract String hgvsProtein();

    @NotNull
    abstract String cDNAPosAndLength();

    @NotNull
    abstract String cdsPosAndLength();

    @NotNull
    public abstract String aaPosAndLength();

    @NotNull
    abstract String distance();

    @NotNull
    abstract String addition();

    // when we use the feature ID it is in practice always a transcript, but this mapping may not hold for every single snpeff annotation!
    public String transcript()
    {
        String transcript = featureID();
        // In case transcripts appear with their version (eg ENST001.1) we strip the version part out.
        if(transcript.contains("."))
        {
            return transcript.substring(0, transcript.indexOf("."));
        }
        else
        {
            return transcript;
        }
    }

    public boolean isTranscriptFeature()
    {
        return featureType().equals(FEATURE_TYPE_TRANSCRIPT);
    }
}
