package com.hartwig.hmftools.apiclients.iclusion.api

import com.hartwig.hmftools.apiclients.iclusion.data.*
import com.hartwig.hmftools.apiclients.iclusion.http.httpClient
import com.squareup.moshi.Moshi
import com.squareup.moshi.kotlin.reflect.KotlinJsonAdapterFactory
import io.reactivex.Observable
import okhttp3.MultipartBody
import okhttp3.OkHttpClient
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import retrofit2.Retrofit
import retrofit2.adapter.rxjava2.RxJava2CallAdapterFactory
import retrofit2.converter.moshi.MoshiConverterFactory

class IclusionApiWrapper(endpoint: String, private val clientId: String, private val clientSecret: String,
                         private val user: String, private val password: String) {
    companion object {
        private fun createApi(httpClient: OkHttpClient, endpoint: String): IclusionApi {
            val retrofit = Retrofit.Builder().baseUrl(endpoint)
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

    private val logger = LogManager.getLogger("IclusionApiWrapper")

    private val httpClient = httpClient()
    private val api = createApi(httpClient, endpoint)
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

    fun studyDetails(): List<IclusionStudyDetails> {
        val studies = studies().blockingIterable().toList()
        logger.info("Queried ${studies.size} studies via iclusion API.")

        logger.info(" Studies with CCMO identifier: ${studies.filterNot { it.ccmo.isEmpty() }.size}")
        logger.info(" Studies with mutations configured: ${studies.filterNot { it.mutations.isEmpty() }.size}")
        logger.info(" Studies with EUDRA identifier: ${studies.filterNot { it.eudra.isEmpty() }.size}")
        logger.info(" Studies with NCT identifier: ${studies.filterNot { it.nct.isEmpty() }.size}")
        logger.info(" Studies with IPN ?: ${studies.filterNot { it.ipn.isNullOrEmpty() }.size}")

        studies.filterNot { it.ipn.isNullOrEmpty() }.forEach { logger.info(it) }

        val indications = indications().blockingIterable().toList().associateBy { it.id }
        logger.info("Queried ${indications.size} indications via iclusion API.")

        val genes = genes().blockingIterable().toList().associateBy { it.id }
        logger.info("Queried ${genes.size} genes via iclusion API.")

        val variants = variants().blockingIterable().toList().associateBy { it.id }
        logger.info("Queried ${variants.size} variants via iclusion API.")

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

    private fun indications(): Observable<IclusionIndication> {
        return api.indications(tokenBearer).flatMapIterable { it }
    }

    private fun genes(): Observable<IclusionGene> {
        return api.genes(tokenBearer).flatMapIterable { it }
    }

    private fun variants(): Observable<IclusionVariant> {
        return api.variants(tokenBearer).flatMapIterable { it }
    }

    private fun studies(): Observable<IclusionStudy> {
        return api.studies(tokenBearer).flatMapIterable { it }
    }
}
