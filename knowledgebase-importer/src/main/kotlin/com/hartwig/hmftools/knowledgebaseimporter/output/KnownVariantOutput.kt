package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.extensions.csv.CsvData

data class KnownVariantOutput(val gene: String, private val transcript: String, private val reference: String,
                              private val chromosome: String, private val position: String, private val ref: String,
                              private val alt: String, private val annotation: String) : CsvData {
    companion object {
        operator fun invoke(transcript: String, reference: String, annotation: String, variant: SomaticVariantEvent): KnownVariantOutput {
            return KnownVariantOutput(
                    variant.gene, transcript, reference, variant.chromosome, variant.position, variant.ref, variant.alt, annotation)
        }
    }
}
