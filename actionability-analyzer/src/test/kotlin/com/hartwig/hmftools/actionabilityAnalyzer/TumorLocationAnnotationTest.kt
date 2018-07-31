package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.FileInputStream

class TumorLocationAnnotationTest : StringSpec() {
    private val logger = LogManager.getLogger("TumorLocationAnnotationTest")
    private val patientDbResourcesPath = "..${File.separator}patient-db${File.separator}src${File.separator}main${File.separator}resources"
    private val tumorLocationMappingCsv = "tumor_location_mapping.csv"

    private val curator = TumorLocationCurator(FileInputStream("$patientDbResourcesPath${File.separator}$tumorLocationMappingCsv"))

    init {
        "all curated primary tumor locations have entries in doid csv" {
            curator.primaryTumorLocations().filter { !ActionabilityAnalyzer.primaryTumorMapping().keys.contains(it.toLowerCase()) }
                    .forEach { logger.error("Missing DOID annotation for $it") }
            curator.primaryTumorLocations().all { ActionabilityAnalyzer.primaryTumorMapping().keys.contains(it.toLowerCase()) } shouldBe true
        }
    }
}
