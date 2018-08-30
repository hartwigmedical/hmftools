package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers.IclusionCDnaAnnotationReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CivicCDnaAnnotationReader : SomaticEventReader<CivicVariantInput, CDnaAnnotation> {
    private const val ENSEMBL_PATTERN = "(ENST[0-9]+)(?:\\.[0-9]+)?"
    private const val ENSEMBL_ANNOTATION_PATTERN = "($ENSEMBL_PATTERN):([^,\\t\\s\\n]+)"

    fun match(event: CivicVariantInput) = ENSEMBL_ANNOTATION_PATTERN.toRegex(RegexOption.IGNORE_CASE).find(event.hgvs_expressions) != null

    override fun read(event: CivicVariantInput): List<CDnaAnnotation> {
        return if (match(event)) {
            return IclusionCDnaAnnotationReader.read(extractAnnotation(event))
        } else {
            emptyList()
        }
    }

    private fun extractAnnotation(event: CivicVariantInput): CivicVariantInput {
        val matchResult = ENSEMBL_ANNOTATION_PATTERN.toRegex(RegexOption.IGNORE_CASE).find(event.hgvs_expressions)!!
        val transcript = matchResult.groupValues[2]
        val annotation = matchResult.groupValues[3]
        return event.copy(representative_transcript = transcript, variant = annotation)
    }
}
