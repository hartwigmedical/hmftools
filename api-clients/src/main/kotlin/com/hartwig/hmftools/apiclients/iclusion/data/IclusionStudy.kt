package com.hartwig.hmftools.apiclients.iclusion.data

data class IclusionStudy(val id: String, val eudra: String, val nct: String, val ipn: String?, val ccmo: String, val title: String,
                         val acronym: String, val indication_ids: List<String>, val mutations: List<IclusionMutation>)
