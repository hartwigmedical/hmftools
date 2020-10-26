package com.hartwig.hmftools.iclusion.api;

import java.util.List;
import java.util.concurrent.TimeUnit;

import com.hartwig.hmftools.iclusion.IclusionCredentials;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
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

public final class IclusionApiMain {

    private static final Logger LOGGER = LogManager.getLogger(IclusionApiMain.class);

    private IclusionApiMain() {
    }

    @NotNull
    public static List<IclusionTrial> readIclusionTrials(@NotNull String iClusionEndpoint, @NotNull IclusionCredentials credentials) {
        OkHttpClient httpClient = buildHttpClient();

        LOGGER.info("Connecting with iClusion API at {} using user '{}'", iClusionEndpoint, credentials.username());
        IclusionApi api = buildIclusionApi(httpClient, iClusionEndpoint);

        LOGGER.info("Requesting iClusion access token");
        String tokenBearer = requestAndBuildTokenBearer(api, credentials);

        LOGGER.info("Querying iClusion trial database using access token");
        List<IclusionObjectStudy> studies = api.studies(tokenBearer).blockingFirst();
        LOGGER.info(" Received {} iClusion studies", studies.size());

        List<IclusionObjectIndication> indications = api.indications(tokenBearer).blockingFirst();
        LOGGER.info(" Received {} iClusion indications", indications.size());

        List<IclusionObjectGene> genes = api.genes(tokenBearer).blockingFirst();
        LOGGER.info(" Received {} iClusion genes", genes.size());

        List<IclusionObjectVariant> variants = api.variants(tokenBearer).blockingFirst();
        LOGGER.info(" Received {} iClusion variants", variants.size());

        LOGGER.info("Closing down connection to iClusion API");
        httpClient.dispatcher().executorService().shutdown();

        List<IclusionTrial> trials = IclusionApiObjectMapper.fromApiObjects(studies, indications, genes, variants);
        LOGGER.info("Created {} trials from iClusion API objects", trials.size());

        return trials;
    }

    @NotNull
    private static OkHttpClient buildHttpClient() {
        Dispatcher dispatcher = new Dispatcher();
        dispatcher.setMaxRequests(100);
        dispatcher.setMaxRequestsPerHost(100);

        return new OkHttpClient.Builder().connectionPool(new ConnectionPool(20, 5, TimeUnit.SECONDS))
                .readTimeout(20, TimeUnit.SECONDS)
                .connectTimeout(20, TimeUnit.SECONDS)
                .writeTimeout(20, TimeUnit.SECONDS)
                .dispatcher(dispatcher)
                .build();
    }

    @NotNull
    private static IclusionApi buildIclusionApi(@NotNull OkHttpClient httpClient, @NotNull String iClusionEndpoint) {
        MoshiConverterFactory moshiConverterFactory =
                MoshiConverterFactory.create(new Moshi.Builder().add(new IclusionResponseAdapter()).build());

        Retrofit retrofit = new Retrofit.Builder().baseUrl(forceTrailingSlash(iClusionEndpoint))
                .addConverterFactory(moshiConverterFactory)
                .addCallAdapterFactory(RxJava2CallAdapterFactory.createAsync())
                .client(httpClient)
                .build();

        return retrofit.create(IclusionApi.class);
    }

    @NotNull
    private static String forceTrailingSlash(@NotNull String iClusionEndpoint) {
        if (!iClusionEndpoint.endsWith("/")) {
            return iClusionEndpoint + "/";
        } else {
            return iClusionEndpoint;
        }
    }

    @NotNull
    private static String requestAndBuildTokenBearer(@NotNull IclusionApi api, @NotNull IclusionCredentials credentials) {
        MultipartBody requestBody = new MultipartBody.Builder().setType(MultipartBody.FORM)
                .addFormDataPart("grant_type", "password")
                .addFormDataPart("client_id", credentials.clientId())
                .addFormDataPart("client_secret", credentials.clientSecret())
                .addFormDataPart("username", credentials.username())
                .addFormDataPart("password", credentials.password())
                .build();

        // "Bearer " has to go in front, see DEV-473
        return "Bearer " + api.requestAccessToken(requestBody).blockingFirst().accessToken;
    }
}
