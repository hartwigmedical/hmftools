package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.common.utils.json.JsonFunctions.string;

import com.google.gson.JsonObject;
import com.hartwig.hmftools.vicc.datamodel.sage.ImmutableSage;
import com.hartwig.hmftools.vicc.datamodel.sage.Sage;

import org.jetbrains.annotations.NotNull;

final class SageObjectFactory {

    private SageObjectFactory() {
    }

    @NotNull
    static Sage create(@NotNull JsonObject sageObject) {
        ViccDatamodelCheckerFactory.sageEntryChecker().check(sageObject);

        return ImmutableSage.builder()
                .gene(string(sageObject, "gene"))
                .entrezId(string(sageObject, "entrez_id"))
                .clinicalManifestation(string(sageObject, "clinical_manifestation"))
                .responseType(string(sageObject, "response_type"))
                .evidenceLabel(string(sageObject, "evidence_label"))
                .drugLabels(string(sageObject, "drug_labels"))
                .germlineOrSomatic(string(sageObject, "germline_or_somatic"))
                .publicationUrl(string(sageObject, "publication_url"))
                .build();
    }
}
