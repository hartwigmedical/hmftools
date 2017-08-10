package com.hartwig.hmftools.civic.api;

import com.hartwig.hmftools.civic.data.CivicGene;
import com.hartwig.hmftools.civic.data.CivicVariant;

import org.jetbrains.annotations.NotNull;

import io.reactivex.Observable;
import retrofit2.http.GET;
import retrofit2.http.Path;

public interface CivicApi {
    @NotNull
    @GET("genes/{id}?identifier_type=entrez_id")
    Observable<CivicGene> getGene(@Path("id") @NotNull final String entrezId);

    @NotNull
    @GET("variants/{id}")
    Observable<CivicVariant> getVariant(@Path("id") @NotNull final String variantId);
}
