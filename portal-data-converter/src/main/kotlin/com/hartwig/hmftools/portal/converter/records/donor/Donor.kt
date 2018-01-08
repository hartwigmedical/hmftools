package com.hartwig.hmftools.portal.converter.records.donor

import com.hartwig.hmftools.portal.converter.extensions.toRecord
import com.hartwig.hmftools.portal.converter.records.Record
import kotlin.reflect.KClass

data class Donor(private val fields: Map<DonorHeader, String>) : Record {
    companion object Factory {
        val header: KClass<DonorHeader> = DonorHeader::class
        private val COUNTRY = "Netherlands"
        operator fun invoke(id: String, gender: String, ageAtEnrollment: String): Donor {
            return Donor(mapOf(
                    DonorHeader.donor_id to id,
                    DonorHeader.donor_sex to gender,
                    DonorHeader.donor_age_at_enrollment to ageAtEnrollment,
                    DonorHeader.donor_region_of_residence to COUNTRY
            ))
        }
    }

    override fun record(): List<String> {
        return fields.toRecord()
    }
}
