package com.hartwig.hmftools.apiclients.diseaseontology.api;

import com.hartwig.hmftools.apiclients.diseaseontology.data.DiseaseOntologyMetadata;

import org.jetbrains.annotations.NotNull;

import io.reactivex.Observable;
import retrofit2.http.GET;
import retrofit2.http.Path;

public interface DiseaseOntologyApi {
    int RETRY_COUNT = 5;

    @NotNull
    @GET("metadata/{id}")
    Observable<DiseaseOntologyMetadata> getMetadataService(@Path("id") final String doid);

    default Observable<DiseaseOntologyMetadata> getMetadata(final String doid) {
        return getMetadataService(doid).retry(RETRY_COUNT);
    }
}
