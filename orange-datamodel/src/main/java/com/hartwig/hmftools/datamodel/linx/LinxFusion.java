package com.hartwig.hmftools.datamodel.linx;

import org.immutables.gson.Gson;
import org.immutables.value.Value;

@Gson.TypeAdapters
@Value.Immutable
public interface LinxFusion {
    String name();
    boolean reported();
    LinxFusionType reportedType();
    FusionPhasedType phased();
    FusionLikelihoodType likelihood();
    int fusedExonUp();
    int fusedExonDown();

    // for orange report
    int chainLinks();
    boolean chainTerminated();
    String domainsKept();
    String domainsLost();

    // for patient report
    String geneStart();
    String geneContextStart();
    String geneTranscriptStart();
    String geneEnd();
    String geneContextEnd();
    String geneTranscriptEnd();
    Double junctionCopyNumber();
}
