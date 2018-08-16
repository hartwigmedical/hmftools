package com.hartwig.hmftools.apiclients.iclusion.data

data class IclusionStudyDetails(val study: IclusionStudy, val indications: List<IclusionIndication>,
                                val mutations: List<IclusionMutationDetails>)
