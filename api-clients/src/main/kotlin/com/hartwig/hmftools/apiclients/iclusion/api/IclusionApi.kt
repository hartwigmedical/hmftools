package com.hartwig.hmftools.apiclients.iclusion.api

import com.hartwig.hmftools.apiclients.iclusion.data.IclusionGene
import com.hartwig.hmftools.apiclients.iclusion.data.IclusionIndication
import com.hartwig.hmftools.apiclients.iclusion.data.IclusionVariant
import com.hartwig.hmftools.apiclients.iclusion.data.Token
import io.reactivex.Observable
import okhttp3.RequestBody
import retrofit2.http.*

interface IclusionApi {
    @POST("oauth/token")
    fun getAccessToken(@Body requestBody: RequestBody): Observable<Token>

    @GET("indications")
    fun indications(@Header("Authorization") tokenBearer: String): Observable<List<IclusionIndication>>

    @GET("indications/{id}")
    fun indication(@Header("Authorization") tokenBearer: String, @Path("id") indicationId: String): Observable<IclusionIndication>

    @GET("genes")
    fun genes(@Header("Authorization") tokenBearer: String): Observable<List<IclusionGene>>

    @GET("variants")
    fun variants(@Header("Authorization") tokenBearer: String): Observable<List<IclusionVariant>>
}
