package com.hartwig.hmftools.knowledgebaseimporter.cgi

class CgiDrug(val drugName: String, drugFamily: String) {
    companion object {
        private fun curate(drugFamily: String): List<String> {
            val family = drugFamily.replace("\\s+".toRegex(), " ")
            return if (family.contains("&")) {
                curateMultiple(family)
            } else {
                listOf(curateSingle(family))
            }
        }

        private fun curateSingle(drugFamily: String): String = when (drugFamily.toLowerCase()) {
            "BCR-ABL inhibitor 1st gen".toLowerCase()                        -> BCR_ABL_INHIBITOR
            "BCR-ABL inhibitor 2nd gen".toLowerCase()                        -> BCR_ABL_INHIBITOR
            "BCR-ABL inhibitor 3rd gen".toLowerCase()                        -> BCR_ABL_INHIBITOR
            "EGFR inhibitor 1st gen".toLowerCase()                           -> EGFR_INHIBITOR
            "EGFR inhibitor 2nd gen".toLowerCase()                           -> EGFR_INHIBITOR
            "EGFR inhibitor 3rd gen".toLowerCase()                           -> EGFR_INHIBITOR
            "3rd generation EGFR inhibitor".toLowerCase()                    -> EGFR_INHIBITOR
            "novel EGFR mAb inhibitor".toLowerCase()                         -> EGFR_MAB_INHIBITOR
            "PD1 Ab".toLowerCase()                                           -> CHECKPOINT_INHIBITOR
            "PD1 Ab inhibitor".toLowerCase()                                 -> CHECKPOINT_INHIBITOR
            "PDL1 inhibitor".toLowerCase()                                   -> CHECKPOINT_INHIBITOR
            "ATP competitive AKT inhibitor".toLowerCase()                    -> AKT_INHIBITOR
            "Allosteric AKT inhibitor".toLowerCase()                         -> AKT_INHIBITOR
            "non-allosteric AKT inhibitor".toLowerCase()                     -> AKT_INHIBITOR
            "Anthracycline antitumor antibiotic".toLowerCase()               -> CHEMOTHERAPY
            "Fluoropyrimidine".toLowerCase()                                 -> CHEMOTHERAPY
            "Taxane".toLowerCase()                                           -> CHEMOTHERAPY
            "Alkylating agent".toLowerCase()                                 -> CHEMOTHERAPY
            "Amylin analogue".toLowerCase()                                  -> CHEMOTHERAPY
            "Purine analog".toLowerCase()                                    -> CHEMOTHERAPY
            "Guanine analog".toLowerCase()                                   -> CHEMOTHERAPY
            "CDK4/CDK6 inhibitor".toLowerCase()                              -> CDK4_6_INHIBITOR
            "PIK3CA inhibitor".toLowerCase()                                 -> PI3K_INHIBITOR
            "PIK3CB inhibitor".toLowerCase()                                 -> PI3K_INHIBITOR
            "PI3K pathway inhibitor".toLowerCase()                           -> PI3K_INHIBITOR
            "PI3K pathway inhibitor (alone or in combination)".toLowerCase() -> PI3K_INHIBITOR
            "novel ALK inhibitor".toLowerCase()                              -> ALK_INHIBITOR
            "JAK inhibitor (alone or in combination)".toLowerCase()          -> JAK_INHIBITOR
            "MEK inhibitor (alone or in combination)".toLowerCase()          -> MEK_INHIBITOR
            "novel MEK inhibitor".toLowerCase()                              -> MEK_INHIBITOR
            "novel FLT3 inhibitor".toLowerCase()                             -> FLT3_INHIBITOR
            "AR inhibitor next gen".toLowerCase()                            -> AR_INHIBITOR
            "novel TRK inhibitor".toLowerCase()                              -> TRK_INHIBITOR
            else                                                             -> drugFamily
        }

        private fun curateMultiple(drugFamily: String): List<String> = when (drugFamily.toLowerCase()) {
            "ERBB2&EGFR inhibitor 2nd gen".toLowerCase() -> listOf(ERBB2_INHIBITOR, EGFR_INHIBITOR)
            "ALK&ROS1 inhibitor".toLowerCase()           -> listOf(ALK_INHIBITOR, ROS1_INHIBITOR)
            else                                         -> drugFamily.split("&").map { curateSingle(it) }
        }
    }

    val drugTypes: List<String> = curate(drugFamily)
}
