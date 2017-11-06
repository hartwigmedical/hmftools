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
    int RETRY_COUNT = 5;

    @NotNull
    @GET("genes/{id}?identifier_type=entrez_id")
    Observable<CivicGene> getGeneService(@Path("id") final int entrezId);

    @NotNull
    @GET("variants/{id}")
    Observable<CivicVariant> getVariantService(@Path("id") final int variantId);

    @NotNull
    @GET("variants")
    Observable<CivicIndexResult<CivicVariant>> getVariantsService(@Query("page") final long page);

    @NotNull
    @GET("evidence_items")
    Observable<CivicIndexResult<CivicEvidenceItem>> getEvidenceItemsService(@Query("page") final long page);

    @NotNull
    @GET("genes")
    Observable<CivicIndexResult<CivicGene>> getGenesService(@Query("page") final long page);

    default Observable<CivicGene> getGene(final int entrezId) {
        return getGeneService(entrezId).retry(RETRY_COUNT);
    }

    default Observable<CivicVariant> getVariant(final int variantId) {
        return getVariantService(variantId).retry(RETRY_COUNT);
    }

    default Observable<CivicIndexResult<CivicVariant>> getVariants(final long page) {
        return getVariantsService(page).retry(RETRY_COUNT);
    }

    default Observable<CivicIndexResult<CivicEvidenceItem>> getEvidenceItems(final long page) {
        return getEvidenceItemsService(page).retry(RETRY_COUNT);
    }

    default Observable<CivicIndexResult<CivicGene>> getGenes(final long page) {
        return getGenesService(page).retry(RETRY_COUNT);
    }
}
