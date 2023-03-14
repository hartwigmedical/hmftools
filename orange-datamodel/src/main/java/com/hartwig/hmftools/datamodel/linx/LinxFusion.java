package com.hartwig.hmftools.datamodel.linx;

import org.immutables.gson.Gson;
import org.immutables.value.Value;

@Gson.TypeAdapters
@Value.Immutable
public abstract class LinxFusion {
    public abstract String name();
    public abstract boolean reported();
    public abstract LinxFusionType reportedType();
    public abstract FusionPhasedType phased();
    public abstract FusionLikelihoodType likelihood();
    public abstract int fusedExonUp();
    public abstract int fusedExonDown();

    // for orange report
    public abstract int chainLinks();
    public abstract boolean chainTerminated();
    public abstract String domainsKept();
    public abstract String domainsLost();

    // for patient report
    public abstract String geneStart();
    public abstract String geneContextStart();
    public abstract String geneTranscriptStart();
    public abstract String geneEnd();
    public abstract String geneContextEnd();
    public abstract String geneTranscriptEnd();
    public abstract Double junctionCopyNumber();
}
