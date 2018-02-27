package com.hartwig.hmftools.apiclients.diseaseontology.api;

import java.util.concurrent.TimeUnit;

import com.google.gson.Gson;
import com.hartwig.hmftools.apiclients.diseaseontology.data.DiseaseOntologyGson;
import com.hartwig.hmftools.apiclients.diseaseontology.data.DiseaseOntologyMetadata;

import org.jetbrains.annotations.NotNull;

import io.reactivex.Observable;
import okhttp3.ConnectionPool;
import okhttp3.Dispatcher;
import okhttp3.OkHttpClient;
import retrofit2.Retrofit;
import retrofit2.adapter.rxjava2.RxJava2CallAdapterFactory;
import retrofit2.converter.gson.GsonConverterFactory;

public class DiseaseOntologyApiWrapper {
    private static final String ENDPOINT = "http://www.disease-ontology.org/api/";
    private final DiseaseOntologyApi api;
    private final OkHttpClient httpClient;
    private final String CANCER_DOID = "DOID:162";

    public DiseaseOntologyApiWrapper() {
        final Dispatcher requestDispatcher = new Dispatcher();
        requestDispatcher.setMaxRequests(100);
        requestDispatcher.setMaxRequestsPerHost(100);
        httpClient = new OkHttpClient.Builder().connectionPool(new ConnectionPool(20, 5, TimeUnit.SECONDS))
                .readTimeout(20, TimeUnit.SECONDS)
                .connectTimeout(20, TimeUnit.SECONDS)
                .writeTimeout(20, TimeUnit.SECONDS)
                .dispatcher(requestDispatcher)
                .build();
        final Gson gson = DiseaseOntologyGson.buildGson();
        final Retrofit retrofit = new Retrofit.Builder().baseUrl(ENDPOINT)
                .addConverterFactory(GsonConverterFactory.create(gson))
                .addCallAdapterFactory(RxJava2CallAdapterFactory.createAsync())
                .client(httpClient)
                .build();
        api = retrofit.create(DiseaseOntologyApi.class);
    }

    @NotNull
    public Observable<DiseaseOntologyMetadata> getMetadata(@NotNull final String doid) {
        return api.getMetadata(doid);
    }

    @NotNull
    public Observable<DiseaseOntologyMetadata> getMetadata(final int doid) {
        return getMetadata("DOID:" + doid);
    }

    @NotNull
    public Observable<String> getAllParentDoids(final String doid) {
        return getMetadata(doid).flatMap(metadata -> Observable.fromIterable(metadata.parentDoids())
                .mergeWith(Observable.fromIterable(metadata.parentDoids())
                        .filter(parentDoid -> !parentDoid.equals(CANCER_DOID))
                        .flatMap(this::getAllParentDoids)));
    }

    @NotNull
    public Observable<String> getAllParentDoids(final int doid) {
        return getAllParentDoids("DOID:" + doid);
    }

    @NotNull
    public Observable<String> getAllChildrenDoids(@NotNull final String doid) {
        return getMetadata(doid).flatMap(metadata -> Observable.fromIterable(metadata.childrenDoids())
                .mergeWith(Observable.fromIterable(metadata.childrenDoids()).flatMap(this::getAllChildrenDoids)));
    }

    @NotNull
    public void releaseResources() {
        httpClient.dispatcher().executorService().shutdown();
    }
}
