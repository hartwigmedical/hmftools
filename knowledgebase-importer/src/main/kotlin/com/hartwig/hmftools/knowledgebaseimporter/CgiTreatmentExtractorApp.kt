package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import com.hartwig.hmftools.extensions.cli.options.filesystem.RequiredInputFileOption
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredOutputOption
import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.extensions.csv.CsvWriter
import com.hartwig.hmftools.knowledgebaseimporter.cgi.CgiActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug
import org.apache.logging.log4j.LogManager
import java.io.File

private val logger = LogManager.getLogger("CgiTreatmentExtractorApp")

fun main(args: Array<String>) {
    val cmd = createOptions().createCommandLine("cgi-treatment-extractor", args)
    val outputFile = File(cmd.getOptionValue(OUTPUT_FILE))
    if (!outputFile.parentFile.exists() && !outputFile.parentFile.mkdirs()) {
        logger.error("Could not create parent directories for $outputFile")
    }
    val cgiActionableRecords = CsvReader.readTSVByName<CgiActionableInput>(cmd.getOptionValue(CGI_BIOMARKERS_LOCATION))
            .mapNotNull { it.corrected() }.map { CgiActionableRecord(it, mapOf()) }
    CsvWriter.writeTSV(bootstrapTreatmentTypeMapping(cgiActionableRecords), "$outputFile")
}

private fun createOptions(): HmfOptions {
    val options = HmfOptions()
    options.add(RequiredInputFileOption(CGI_BIOMARKERS_LOCATION, "path to cgi biomarkers file"))
    options.add(RequiredOutputOption(OUTPUT_FILE, "path to output file"))
    return options
}

private fun bootstrapTreatmentTypeMapping(records: List<CgiActionableRecord>): List<HmfDrug> {
    return records.flatMap { it.cgiDrugs }.groupBy { it.drugName }
            .mapValues { (key, value) ->
                val drugFamilies = value.flatMap { it.drugTypes }.toSet()
                if (drugFamilies.size > 1) drugFamilies.filterNot { it == key } else drugFamilies
            }.map { HmfDrug(it.key, it.value.joinToString(";")) }
}
