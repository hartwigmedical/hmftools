package com.hartwig.hmftools.portal.converter.records.specimen

import com.hartwig.hmftools.portal.converter.extensions.toRecord
import com.hartwig.hmftools.portal.converter.records.Record
import kotlin.reflect.KClass

data class Specimen(private val fields: Map<SpecimenHeader, String>) : Record {
    companion object Factory {
        val header: KClass<SpecimenHeader> = SpecimenHeader::class
        private val STORAGE = "1"
        operator fun invoke(id: String, type: String, typeOther: String, donorId: String): Specimen {
            return Specimen(mapOf(
                    SpecimenHeader.specimen_id to id,
                    SpecimenHeader.specimen_type to type,
                    SpecimenHeader.specimen_type_other to typeOther,
                    SpecimenHeader.donor_id to donorId,
                    SpecimenHeader.specimen_storage to STORAGE
            ))
        }
    }

    override fun record(): List<String> {
        return fields.toRecord()
    }
}
