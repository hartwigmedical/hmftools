package com.hartwig.hmftools.apiclients.iclusion.api

import com.hartwig.hmftools.apiclients.iclusion.data.*
import com.hartwig.hmftools.apiclients.iclusion.http.httpClient
import com.squareup.moshi.Moshi
import com.squareup.moshi.kotlin.reflect.KotlinJsonAdapterFactory
import io.reactivex.Observable
import okhttp3.MultipartBody
import okhttp3.OkHttpClient
import retrofit2.Retrofit
import retrofit2.adapter.rxjava2.RxJava2CallAdapterFactory
import retrofit2.converter.moshi.MoshiConverterFactory

private const val ICLUSION_ENDPOINT = "http://api.iclusion.nl/"

class IclusionApiWrapper {
    companion object {
        private fun createApi(httpClient: OkHttpClient): IclusionApi {
            val retrofit = Retrofit.Builder().baseUrl(ICLUSION_ENDPOINT)
                    .addConverterFactory(moshiConverter())
                    .addCallAdapterFactory(RxJava2CallAdapterFactory.createAsync())
                    .client(httpClient)
                    .build()
            return retrofit.create(IclusionApi::class.java)
        }

        private fun moshiConverter(): MoshiConverterFactory {
            return MoshiConverterFactory.create(Moshi.Builder().add(KotlinJsonAdapterFactory()).add(IclusionResponseAdapter).build())
        }
    }

    private val httpClient = httpClient()
    private val api = createApi(httpClient)

    fun getAccessToken(clientId: String, clientSecret: String, user: String, password: String): Observable<Token> {
        val requestBody = MultipartBody.Builder().setType(MultipartBody.FORM)
                .addFormDataPart("grant_type", "password")
                .addFormDataPart("client_id", clientId)
                .addFormDataPart("client_secret", clientSecret)
                .addFormDataPart("username", user)
                .addFormDataPart("password", password)
                .build()
        return api.getAccessToken(requestBody)
    }

    fun indications(token: Token): Observable<IclusionIndication> {
        return api.indications(tokenBearer(token)).flatMapIterable { it }
    }

    fun indication(token: Token, indicationId: String): Observable<IclusionIndication> {
        return api.indication(tokenBearer(token), indicationId)
    }

    fun genes(token: Token): Observable<IclusionGene> {
        return api.genes(tokenBearer(token)).flatMapIterable { it }
    }

    fun variants(token: Token): Observable<IclusionVariant> {
        return api.variants(tokenBearer(token)).flatMapIterable { it }
    }

    fun studies(token: Token): Observable<IclusionStudy> {
        return api.studies(tokenBearer(token)).flatMapIterable { it }
    }

    fun studyDetails(token: Token): List<IclusionStudyDetails> {
        val studies = studies(token).blockingIterable().toList()
        val indications = indications(token).blockingIterable().toList().associateBy { it.id }
        val genes = genes(token).blockingIterable().toList().associateBy { it.id }
        val variants = variants(token).blockingIterable().toList().associateBy { it.id }
        return studies.filterNot { it.mutations.isEmpty() }.map {
            val indicationsForStudy = it.indication_ids.mapNotNull { indications[it] }
            val mutationsDetails = it.mutations.mapNotNull {
                val geneName = genes[it.gene_id]?.gene_name
                val variantName = variants[it.variant_id]?.variant_name
                geneName ?: return@mapNotNull null
                variantName ?: return@mapNotNull null
                return@mapNotNull IclusionMutationDetails(geneName, variantName)
            }
            IclusionStudyDetails(it, indicationsForStudy, mutationsDetails)
        }.filterNot { it.mutations.isEmpty() }
    }

    fun close() {
        httpClient.dispatcher().executorService().shutdown()
    }

    private fun tokenBearer(token: Token): String = "Bearer ${token.access_token}"
}
