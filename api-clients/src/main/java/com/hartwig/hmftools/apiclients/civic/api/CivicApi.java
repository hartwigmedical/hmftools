package com.hartwig.hmftools.apiclients.civic.api;

import com.hartwig.hmftools.apiclients.civic.data.CivicEvidenceItem;
import com.hartwig.hmftools.apiclients.civic.data.CivicGene;
import com.hartwig.hmftools.apiclients.civic.data.CivicIndexResult;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;

import org.jetbrains.annotations.NotNull;

import io.reactivex.Observable;
import retrofit2.http.GET;
import retrofit2.http.Path;
import retrofit2.http.Query;

public interface CivicApi {
    @NotNull
    @GET("genes/{id}?identifier_type=entrez_id")
    Observable<CivicGene> getGene(@Path("id") final int entrezId);

    @NotNull
    @GET("variants/{id}")
    Observable<CivicVariant> getVariant(@Path("id") final int variantId);

    @NotNull
    @GET("variants")
    Observable<CivicIndexResult<CivicVariant>> getVariants(@Query("page") final long page);

    @NotNull
    @GET("evidence_items")
    Observable<CivicIndexResult<CivicEvidenceItem>> getEvidenceItems(@Query("page") final long page);

    @NotNull
    @GET("genes")
    Observable<CivicIndexResult<CivicGene>> getGenes(@Query("page") final long page);
}
