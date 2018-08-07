package com.hartwig.hmftools.apiclients.iclusion.data

import com.squareup.moshi.FromJson

object IclusionResponseAdapter {
    @FromJson
    fun indications(json: Map<String, Indication>): List<Indication> {
        return json.values.toList()
    }
}
