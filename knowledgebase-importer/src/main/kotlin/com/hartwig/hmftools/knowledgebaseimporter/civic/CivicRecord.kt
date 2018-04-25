package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.google.common.collect.Multimap
import org.apache.commons.csv.CSVRecord
import java.util.regex.Pattern

class CivicRecord(variantCSVRecord: CSVRecord, variantEvidenceMap: Multimap<String, CivicEvidence>) {
    val transcript: String = variantCSVRecord["representative_transcript"] ?: "-"
    val variant: String = variantCSVRecord["variant"]
    val gene: String = variantCSVRecord["gene"]
    val chromosome: String = variantCSVRecord["chromosome"].orEmpty()
    val start: String = variantCSVRecord["start"].orEmpty()
    val stop: String = variantCSVRecord["stop"].orEmpty()
    val ref: String = variantCSVRecord["reference_bases"].orEmpty()
    val alt: String = variantCSVRecord["variant_bases"].orEmpty()
    val evidence: Collection<CivicEvidence> = variantEvidenceMap[variantCSVRecord["variant_id"]]
    val hgvs: String = extractHgvs(variantCSVRecord["hgvs_expressions"].orEmpty())

    private fun extractHgvs(hgvsField: String): String {
        val matcher = Pattern.compile("(ENST[0-9]+\\.[0-9+]:c\\.[0-9][^,\\t\\s\\n]+)").matcher(hgvsField)
        return if (matcher.find()) {
            matcher.group(1)
        } else {
            ""
        }
    }
}
