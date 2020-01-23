package com.hartwig.hmftools.iclusion.api;

import java.util.List;
import java.util.concurrent.TimeUnit;

import com.squareup.moshi.Moshi;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import okhttp3.ConnectionPool;
import okhttp3.Dispatcher;
import okhttp3.MultipartBody;
import okhttp3.OkHttpClient;
import retrofit2.Retrofit;
import retrofit2.adapter.rxjava2.RxJava2CallAdapterFactory;
import retrofit2.converter.moshi.MoshiConverterFactory;

public final class IclusionApiWrapper {

    private static final Logger LOGGER = LogManager.getLogger(IclusionApiWrapper.class);

    private IclusionApiWrapper() {
    }

    public static void connectWithIclusionApi(@NotNull String iClusionEndpoint, @NotNull String iClusionClientId,
            @NotNull String iClusionClientSecret, @NotNull String iClusionUsername, @NotNull String iClusionPassword) {
        LOGGER.info("Connecting with iClusion API on {}", iClusionEndpoint);

        Dispatcher dispatcher = new Dispatcher();
        dispatcher.setMaxRequests(100);
        dispatcher.setMaxRequestsPerHost(100);

        OkHttpClient client = new OkHttpClient.Builder().connectionPool(new ConnectionPool(20, 5, TimeUnit.SECONDS))
                .readTimeout(20, TimeUnit.SECONDS)
                .connectTimeout(20, TimeUnit.SECONDS)
                .writeTimeout(20, TimeUnit.SECONDS)
                .dispatcher(dispatcher)
                .build();

        MoshiConverterFactory moshiConverterFactory =
                MoshiConverterFactory.create(new Moshi.Builder().add(new IclusionResponseAdapter()).build());

        Retrofit retrofit = new Retrofit.Builder().baseUrl(iClusionEndpoint)
                .addConverterFactory(moshiConverterFactory)
                .addCallAdapterFactory(RxJava2CallAdapterFactory.createAsync())
                .client(client)
                .build();

        IclusionApi api = retrofit.create(IclusionApi.class);

        MultipartBody requestBody = new MultipartBody.Builder().setType(MultipartBody.FORM)
                .addFormDataPart("grant_type", "password")
                .addFormDataPart("client_id", iClusionClientId)
                .addFormDataPart("client_secret", iClusionClientSecret)
                .addFormDataPart("username", iClusionUsername)
                .addFormDataPart("password", iClusionPassword)
                .build();

        String tokenBearer = "Bearer " + api.requestAccessToken(requestBody).blockingFirst().accessToken;
        System.out.println("TokenBearer = " + tokenBearer);

        List<IclusionStudy> studies = api.studies(tokenBearer).blockingFirst();
        System.out.println("Studies = " + studies.size());

        List<IclusionIndication> indications = api.indications(tokenBearer).blockingFirst();
        System.out.println("Indications = " + indications.size());

        List<IclusionGene> genes = api.genes(tokenBearer).blockingFirst();
        System.out.println("Genes = " + genes.size());

        List<IclusionVariant> variants = api.variants(tokenBearer).blockingFirst();
        System.out.println("Variants = " + variants.size());

        client.dispatcher().executorService().shutdown();
    }
}
