package com.hartwig.hmftools.portal.converter.extensions

import com.hartwig.hmftools.common.ecrf.projections.PortalClinicalData
import com.hartwig.hmftools.portal.converter.SampleClinicalData
import org.apache.logging.log4j.LogManager
import java.io.File

private val logger = LogManager.getLogger("FileExtensions")

fun File.readSamplesData(): List<SampleClinicalData> {
    return PortalClinicalData.readRecords(this.canonicalPath).map { patientData ->
        val sampleData = SampleClinicalData(patientData)
        if (sampleData.cpctId.isBlank()) {
            logger.warn("Read a CSV record with null CPCT ID")
        }
        sampleData
    }
}
