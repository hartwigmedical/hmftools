package com.hartwig.hmftools.portal.converter.records.sample

import com.hartwig.hmftools.portal.converter.extensions.toRecord
import com.hartwig.hmftools.portal.converter.records.Record
import kotlin.reflect.KClass

data class Sample(private val fields: Map<SampleHeader, String>) : Record {
    companion object Factory {
        val header: KClass<SampleHeader> = SampleHeader::class
        operator fun invoke(id: String): Sample {
            return Sample(mapOf(
                    SampleHeader.analyzed_sample_id to id,
                    SampleHeader.specimen_id to id
            ))
        }
    }

    override fun record(): List<String> {
        return fields.toRecord()
    }
}
