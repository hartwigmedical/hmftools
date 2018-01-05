package com.hartwig.hmftools.portal.converter.donor

import com.hartwig.hmftools.portal.converter.DEFAULT_VALUE
import com.hartwig.hmftools.portal.converter.Record
import com.hartwig.hmftools.portal.converter.extensions.toRecord
import kotlin.reflect.KClass

data class Donor(private val fields: Map<DonorHeader, String>) : Record {
    companion object Factory {
        val header: KClass<DonorHeader> = DonorHeader::class
        private val COUNTRY = "Netherlands"
        operator fun invoke(id: String, gender: Int, ageAtEnrollment: Int?): Donor {
            return Donor(mapOf(
                    DonorHeader.donor_id to id,
                    DonorHeader.donor_sex to gender.toString(),
                    DonorHeader.donor_age_at_enrollment to (ageAtEnrollment?.toString() ?: DEFAULT_VALUE),
                    DonorHeader.donor_region_of_residence to COUNTRY
            ))
        }
    }

    override fun record(): List<String> {
        return fields.toRecord()
    }
}
