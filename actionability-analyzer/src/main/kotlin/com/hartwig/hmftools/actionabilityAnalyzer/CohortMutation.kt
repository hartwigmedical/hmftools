package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.extensions.csv.CsvData

data class CohortMutation(val sampleId: String, val chromosome: String, val position: String, val ref: String, val alt: String,
                          val type: String, val gene: String, val impact: String, val worstCodingEffect: String,
                          val canonicalCodingEffect: String, val transcriptId: String) : CsvData
