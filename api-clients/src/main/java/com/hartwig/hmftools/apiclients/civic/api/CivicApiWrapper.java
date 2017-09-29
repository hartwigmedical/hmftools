package com.hartwig.hmftools.apiclients.civic.api;

import java.util.concurrent.TimeUnit;
import java.util.function.Function;

import com.google.gson.Gson;
import com.hartwig.hmftools.apiclients.civic.data.CivicApiDataGson;
import com.hartwig.hmftools.apiclients.civic.data.CivicApiMetadata;
import com.hartwig.hmftools.apiclients.civic.data.CivicEvidenceItem;
import com.hartwig.hmftools.apiclients.civic.data.CivicGene;
import com.hartwig.hmftools.apiclients.civic.data.CivicIndexResult;
import com.hartwig.hmftools.apiclients.civic.data.CivicVariant;
import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;

import io.reactivex.Observable;
import okhttp3.ConnectionPool;
import okhttp3.Dispatcher;
import okhttp3.OkHttpClient;
import retrofit2.Retrofit;
import retrofit2.adapter.rxjava2.RxJava2CallAdapterFactory;
import retrofit2.converter.gson.GsonConverterFactory;

public class CivicApiWrapper {
    private static final String CIVIC_API_ENDPOINT = "https://civic.genome.wustl.edu/api/";
    private final CivicApi api;
    private final OkHttpClient httpClient;

    public CivicApiWrapper() {
        final Dispatcher requestDispatcher = new Dispatcher();
        requestDispatcher.setMaxRequests(100);
        requestDispatcher.setMaxRequestsPerHost(100);
        httpClient = new OkHttpClient.Builder().connectionPool(new ConnectionPool(20, 5, TimeUnit.SECONDS))
                .readTimeout(20, TimeUnit.SECONDS)
                .dispatcher(requestDispatcher)
                .build();
        final Gson gson = CivicApiDataGson.buildGson();
        final Retrofit retrofit = new Retrofit.Builder().baseUrl(CIVIC_API_ENDPOINT)
                .addConverterFactory(GsonConverterFactory.create(gson))
                .addCallAdapterFactory(RxJava2CallAdapterFactory.createAsync())
                .client(httpClient)
                .build();
        api = retrofit.create(CivicApi.class);
    }

    @NotNull
    public Observable<CivicVariant> getVariantsForGene(final int entrezId) {
        return api.getGene(entrezId)
                .flatMapIterable(CivicGene::variantMetadatas)
                .flatMap(variantMetadata -> api.getVariant(variantMetadata.id()));
    }

    @NotNull
    public Observable<CivicVariant> getVariantsContaining(final int entrezId, @NotNull final Variant variant) {
        return getVariantsForGene(entrezId).filter(civicVariant -> civicVariant.coordinates().containsVariant(variant));
    }

    @NotNull
    public Observable<CivicVariant> getVariantMatches(final int entrezId, @NotNull final Variant variant) {
        return getVariantsForGene(entrezId).filter(civicVariant -> civicVariant.coordinates().equals(variant));
    }

    @NotNull
    public Observable<CivicVariant> getAllVariants() {
        return getAllFromPaginatedEndpoint(api::getVariants).flatMap(variant -> api.getVariant(variant.id()));
    }

    @NotNull
    public Observable<CivicEvidenceItem> getAllEvidenceItems() {
        return getAllFromPaginatedEndpoint(api::getEvidenceItems);
    }

    @NotNull
    public Observable<CivicGene> getAllGenes() {
        return getAllFromPaginatedEndpoint(api::getGenes);
    }

    @NotNull
    private <T> Observable<T> getAllFromPaginatedEndpoint(@NotNull final Function<Long, Observable<CivicIndexResult<T>>> endpoint) {
        return endpoint.apply(1L).flatMap(indexResult -> {
            final CivicApiMetadata metadata = indexResult.meta();
            final Observable<T> firstPageResults = Observable.fromIterable(indexResult.records());
            final Observable<T> nextPagesResults = Observable.range(2, metadata.totalPages() - 1)
                    .flatMap(page -> endpoint.apply((long) page).flatMapIterable(CivicIndexResult::records));
            return firstPageResults.mergeWith(nextPagesResults);
        });
    }

    public void releaseResources() {
        httpClient.dispatcher().executorService().shutdown();
    }
}
