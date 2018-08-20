package com.hartwig.hmftools.apiclients.iclusion.data

import com.squareup.moshi.FromJson

object IclusionResponseAdapter {
    @FromJson
    fun indications(json: Map<String, IclusionIndication>): List<IclusionIndication> {
        return json.values.toList()
    }

    @FromJson
    fun genes(json: Map<String, IclusionGene>): List<IclusionGene> {
        return json.values.toList()
    }

    @FromJson
    fun variants(json: Map<String, IclusionVariant>): List<IclusionVariant> {
        return json.values.toList()
    }

    @FromJson
    fun studies(json: Map<String, IclusionStudy>): List<IclusionStudy> {
        return json.values.toList()
    }
}
