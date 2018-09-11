package com.hartwig.hmftools.apiclients.iclusion.data

data class IclusionIndication(val id: String, val doid: String?, val doid2: String?, val indication_name: String,
                              val indication_name_full: String, val parent_id: String, val node_ids: List<String>) {

    val doidSet: Set<String> = listOfNotNull(doid, doid2).toSet()
}
