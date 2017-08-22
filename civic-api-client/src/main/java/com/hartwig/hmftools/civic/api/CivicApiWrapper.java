package com.hartwig.hmftools.civic.api;

import com.google.gson.Gson;
import com.hartwig.hmftools.civic.data.CivicApiDataGson;
import com.hartwig.hmftools.civic.data.CivicGene;
import com.hartwig.hmftools.civic.data.CivicVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

import io.reactivex.Observable;
import okhttp3.OkHttpClient;
import retrofit2.Retrofit;
import retrofit2.adapter.rxjava2.RxJava2CallAdapterFactory;
import retrofit2.converter.gson.GsonConverterFactory;

public class CivicApiWrapper {
    private static final String CIVIC_API_ENDPOINT = "https://civic.genome.wustl.edu/api/";

    private static final CivicApi civicApi;

    static {
        OkHttpClient client = new OkHttpClient.Builder().build();
        final Gson gson = CivicApiDataGson.buildGson();
        final Retrofit retrofit = new Retrofit.Builder().baseUrl(CIVIC_API_ENDPOINT)
                .addConverterFactory(GsonConverterFactory.create(gson))
                .addCallAdapterFactory(RxJava2CallAdapterFactory.create())
                .client(client)
                .build();
        civicApi = retrofit.create(CivicApi.class);
    }

    public static Observable<CivicGene> getGene(final int entrezId) {
        return civicApi.getGene(entrezId);
    }

    public static Observable<CivicVariant> getVariantsForGene(final int entrezId) {
        return civicApi.getGene(entrezId)
                .flatMapIterable(CivicGene::variantMetadatas)
                .filter(variantMetadata -> variantMetadata.acceptedEvidenceItems() > 0)
                .flatMap(variantMetadata -> civicApi.getVariant(variantMetadata.id()));
    }

    public static Observable<CivicVariant> getVariantsAtPosition(final int entrezId, @NotNull final SomaticVariant somaticVariant) {
        return getVariantsForGene(entrezId).filter(variant -> variant.coordinates().containsVariant(somaticVariant));
    }

    public static Observable<CivicVariant> getVariantMatches(final int entrezId, @NotNull final SomaticVariant somaticVariant) {
        return getVariantsForGene(entrezId).filter(variant -> variant.coordinates().equals(somaticVariant));
    }
}
