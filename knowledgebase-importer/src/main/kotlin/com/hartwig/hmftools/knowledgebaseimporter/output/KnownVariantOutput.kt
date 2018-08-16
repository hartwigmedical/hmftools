package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.extensions.csv.CsvData

data class KnownVariantOutput(val gene: String, private val transcript: String, private val info: String,
                              private val chromosome: String, private val position: String, private val ref: String,
                              private val alt: String) : CsvData {
    companion object {
        operator fun invoke(transcript: String, info: String, variant: SomaticVariantEvent): KnownVariantOutput {
            return KnownVariantOutput(variant.gene, transcript, info, variant.chromosome, variant.position, variant.ref, variant.alt)
        }
    }
}
