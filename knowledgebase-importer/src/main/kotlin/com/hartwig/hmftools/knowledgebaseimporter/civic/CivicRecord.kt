package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.google.common.collect.Multimap
import org.apache.commons.csv.CSVRecord
import java.util.regex.Pattern

data class CivicRecord(val transcript: String, val gene: String, val chromosome: String, val start: String, val stop: String,
                       val ref: String, val alt: String, val evidence: Collection<CivicEvidence>, val hgvs: String) {
    companion object Factory {
        private val hgvsCDnaPattern = Pattern.compile("(ENST[0-9]+\\.[0-9+]:c\\.[0-9][^,\\t\\s\\n]+)")

        operator fun invoke(variantCSVRecord: CSVRecord, variantEvidenceMap: Multimap<String, CivicEvidence>): CivicRecord {
            val hgvs = extractHgvs(variantCSVRecord.get("hgvs_expressions"))
            val evidence = variantEvidenceMap[variantCSVRecord.get("variant_id")]
            val representativeTranscript = variantCSVRecord.get("representative_transcript")
            val transcript = if (representativeTranscript.isEmpty()) {
                "-"
            } else {
                representativeTranscript
            }
            return CivicRecord(transcript,
                               variantCSVRecord.get("gene"),
                               variantCSVRecord.get("chromosome"),
                               variantCSVRecord.get("start"),
                               variantCSVRecord.get("stop"),
                               variantCSVRecord.get("reference_bases"),
                               variantCSVRecord.get("variant_bases"),
                               evidence,
                               hgvs)
        }

        private fun extractHgvs(hgvsField: String): String {
            val matcher = hgvsCDnaPattern.matcher(hgvsField)
            return if (matcher.find()) {
                matcher.group(1)
            } else {
                ""
            }
        }
    }
}
