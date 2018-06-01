package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.extensions.csv.CsvData

data class KnownVariantOutput(private val transcript: String, private val additionalInfo: String,
                              private val variant: SomaticVariantEvent) : CsvData