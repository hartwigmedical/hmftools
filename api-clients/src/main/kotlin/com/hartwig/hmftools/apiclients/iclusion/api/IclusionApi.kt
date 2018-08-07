package com.hartwig.hmftools.apiclients.iclusion.api

import com.hartwig.hmftools.apiclients.iclusion.data.Indication
import com.hartwig.hmftools.apiclients.iclusion.data.Token
import io.reactivex.Observable
import okhttp3.RequestBody
import retrofit2.http.*

interface IclusionApi {
    @Headers("Accept: application/json")
    @POST("oauth/token")
    fun getAccessToken(@Body requestBody: RequestBody): Observable<Token>

    @Headers("Accept: application/json")
    @GET("indications")
    fun indications(@Header("Authorization") tokenBearer: String): Observable<List<Indication>>
}
