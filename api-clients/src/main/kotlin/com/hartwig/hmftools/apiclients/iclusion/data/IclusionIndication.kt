package com.hartwig.hmftools.apiclients.iclusion.data

data class IclusionIndication(val id: String, val doid: String?, val doid2: String?, val indication_name: String,
                              val indication_name_full: String,
                              val node_ids: List<String>)
