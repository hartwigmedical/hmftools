package com.hartwig.hmftools.iclusion.api;

import java.util.List;

import org.jetbrains.annotations.NotNull;

import io.reactivex.Observable;
import okhttp3.RequestBody;
import retrofit2.http.Body;
import retrofit2.http.GET;
import retrofit2.http.Header;
import retrofit2.http.POST;

public interface IclusionApi {

    @POST("oauth/token")
    @NotNull
    Observable<IclusionToken> requestAccessToken(@Body RequestBody requestBody);

    @GET("indications")
    @NotNull
    Observable<List<IclusionObjectIndication>> indications(@Header("Authorization") @NotNull String tokenBearer);

    @GET("genes")
    @NotNull
    Observable<List<IclusionObjectGene>> genes(@Header("Authorization") @NotNull String tokenBearer);

    @GET("variants")
    @NotNull
    Observable<List<IclusionObjectVariant>> variants(@Header("Authorization") @NotNull String tokenBearer);

    @GET("studies")
    @NotNull
    Observable<List<IclusionObjectStudy>> studies(@Header("Authorization") @NotNull String tokenBearer);
}
