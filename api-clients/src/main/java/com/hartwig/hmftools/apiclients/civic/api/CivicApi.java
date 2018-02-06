package com.hartwig.hmftools.apiclients.civic.api;

import com.hartwig.hmftools.apiclients.civic.data.CivicEvidenceItem;
import com.hartwig.hmftools.apiclients.civic.data.CivicGene;
import com.hartwig.hmftools.apiclients.civic.data.CivicIndexResult;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariantData;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariantWithEvidence;

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
    Observable<CivicVariantWithEvidence> getVariantService(@Path("id") final int variantId);

    @NotNull
    @GET("variants")
    Observable<CivicIndexResult<CivicVariantData>> getVariantsService(@Query("page") final long page, @Query("count") final long count);

    @NotNull
    @GET("evidence_items")
    Observable<CivicIndexResult<CivicEvidenceItem>> getEvidenceItemsService(@Query("page") final long page,
            @Query("count") final long count);

    @NotNull
    @GET("genes")
    Observable<CivicIndexResult<CivicGene>> getGenesService(@Query("page") final long page, @Query("count") final long count);

    default Observable<CivicGene> getGene(final int entrezId) {
        return getGeneService(entrezId).retry(RETRY_COUNT);
    }

    default Observable<CivicVariantWithEvidence> getVariant(final int variantId) {
        return getVariantService(variantId).retry(RETRY_COUNT);
    }

    default Observable<CivicIndexResult<CivicVariantData>> getVariants(final long page, final long count) {
        return getVariantsService(page, count).retry(RETRY_COUNT);
    }

    default Observable<CivicIndexResult<CivicEvidenceItem>> getEvidenceItems(final long page, final long count) {
        return getEvidenceItemsService(page, count).retry(RETRY_COUNT);
    }

    default Observable<CivicIndexResult<CivicGene>> getGenes(final long page, final long count) {
        return getGenesService(page, count).retry(RETRY_COUNT);
    }
}
