package com.hartwig.hmftools.knowledgebaseimporter.cgi

import org.apache.commons.csv.CSVRecord
import org.apache.logging.log4j.LogManager

private val logger = LogManager.getLogger("CgiDrugReader")

fun readCgiDrugs(record: CSVRecord): List<CgiDrug> {
    val familyField = record["Drug family"].orEmpty()
    val namesField = record["Drug"].orEmpty()
    val cgiDrugs = readGenericDrugs(familyField, namesField) + readGroupedDrugs(familyField, namesField) +
            readDrugs(familyField, namesField)
    if (cgiDrugs.isEmpty()) {
        logger.error("Encountered unexpected formats when reading cgi drugs and drug families.")
        logger.error("Drug names field: $namesField")
        logger.error("Drug family field: $familyField")
    }
    return cgiDrugs
}

private fun readGenericDrugs(familyField: String, namesField: String): List<CgiDrug> {
    return if (!hasContent(namesField)) {
        val drugFamilies = familyField.removeSurrounding("[", "]").split(";")
        return drugFamilies.map { CgiDrug(it, it) }
    } else {
        listOf()
    }
}

private fun readGroupedDrugs(familyField: String, namesField: String): List<CgiDrug> {
    return if (isGroupedDrugFamily(familyField) && hasContent(namesField) && namesField.split(",").isNotEmpty()) {
        val drugFamily = familyField.removeSurrounding("[", "]")
        val drugNames = namesField.removeSurrounding("[", "]").split(",")
                .map { it.replace(";", " + ") }.filterNot { it.isBlank() }
        drugNames.map { CgiDrug(it, drugFamily) }
    } else {
        listOf()
    }
}

private fun readDrugs(familyField: String, namesField: String): List<CgiDrug> {
    return if (isNotGroupedFamily(familyField) && hasContent(namesField)) {
        val drugFamilies = familyField.split(";")
        val drugNames = namesField.split(";")
        if (drugFamilies.size != drugNames.size) {
            logger.error("Drug name and family fields did not have the same number of elements after split")
            logger.error("Drug names field: $namesField")
            logger.error("Drug family field: $familyField")
            listOf()
        } else {
            drugNames.zip(drugFamilies).map { CgiDrug(it.first, it.second) }
        }
    } else {
        listOf()
    }
}

private fun hasContent(field: String): Boolean {
    return field.removeSurrounding("[", "]").isNotBlank()
}

private fun isGroupedDrugFamily(drugFamilyField: String): Boolean {
    return drugFamilyField.startsWith("[") && drugFamilyField.endsWith("]") && !drugFamilyField.contains(",")
            && !drugFamilyField.contains(";")
}

private fun isNotGroupedFamily(drugFamilyField: String): Boolean {
    return !drugFamilyField.startsWith("[") && !drugFamilyField.endsWith("]") && !drugFamilyField.contains(",")
}
