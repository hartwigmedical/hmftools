package com.hartwig.hmftools.apiclients.diseaseontology.api;

import com.hartwig.hmftools.apiclients.diseaseontology.data.DiseaseOntologyMetadata;

import org.jetbrains.annotations.NotNull;

import io.reactivex.Observable;
import retrofit2.http.GET;
import retrofit2.http.Path;

public interface DiseaseOntologyApi {
    @NotNull
    @GET("metadata/{id}")
    Observable<DiseaseOntologyMetadata> getMetadata(@Path("id") final String doid);
}
