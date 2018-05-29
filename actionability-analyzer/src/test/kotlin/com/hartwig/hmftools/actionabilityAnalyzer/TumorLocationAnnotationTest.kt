package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec
import java.io.File
import java.io.FileInputStream


class TumorLocationAnnotationTest : StringSpec() {
    private val patientDbResourcesPath = "..${File.separator}patient-db${File.separator}src${File.separator}main${File.separator}resources"
    private val tumorLocationMappingCsv = "tumor_location_mapping.csv"

    private val curator = TumorLocationCurator(FileInputStream("$patientDbResourcesPath${File.separator}$tumorLocationMappingCsv"))

    init {
        "curated primary tumor locations match annotated primary tumor locations" {
            curator.primaryTumorLocations() shouldBe ActionabilityAnalyzer.primaryTumorMapping().keys
        }

    }
}
