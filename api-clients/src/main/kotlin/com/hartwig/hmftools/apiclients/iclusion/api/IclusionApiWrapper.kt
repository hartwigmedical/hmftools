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

class IclusionApiWrapper(private val clientId: String, private val clientSecret: String, private val user: String,
                         private val password: String) {
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
    private val tokenBearer = "Bearer ${getAccessToken().blockingFirst().access_token}"

    private fun getAccessToken(): Observable<Token> {
        val requestBody = MultipartBody.Builder().setType(MultipartBody.FORM)
                .addFormDataPart("grant_type", "password")
                .addFormDataPart("client_id", clientId)
                .addFormDataPart("client_secret", clientSecret)
                .addFormDataPart("username", user)
                .addFormDataPart("password", password)
                .build()
        return api.getAccessToken(requestBody)
    }

    fun indications(): Observable<IclusionIndication> {
        return api.indications(tokenBearer).flatMapIterable { it }
    }

    fun indication(indicationId: String): Observable<IclusionIndication> {
        return api.indication(tokenBearer, indicationId)
    }

    fun genes(): Observable<IclusionGene> {
        return api.genes(tokenBearer).flatMapIterable { it }
    }

    fun variants(): Observable<IclusionVariant> {
        return api.variants(tokenBearer).flatMapIterable { it }
    }

    fun studies(): Observable<IclusionStudy> {
        return api.studies(tokenBearer).flatMapIterable { it }
    }

    fun studyDetails(): List<IclusionStudyDetails> {
        val studies = studies().blockingIterable().toList()
        val indications = indications().blockingIterable().toList().associateBy { it.id }
        val genes = genes().blockingIterable().toList().associateBy { it.id }
        val variants = variants().blockingIterable().toList().associateBy { it.id }
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
}
