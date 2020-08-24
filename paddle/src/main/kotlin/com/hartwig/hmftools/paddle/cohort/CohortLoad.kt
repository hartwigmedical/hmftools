package com.hartwig.hmftools.paddle.cohort

import com.hartwig.hmftools.common.drivercatalog.dnds.DndsMutationalLoadFile


data class CohortLoad(val cohortSize: Int, val indel: Int, val snv: Int) {

    companion object {
        fun fromFile(file: String): CohortLoad {
            val dndsSampleLoad = DndsMutationalLoadFile.read(file)
            var indelLoad = 0
            var snvLoad = 0
            for (dndsMutationalLoad in dndsSampleLoad) {
                indelLoad += dndsMutationalLoad.indelTotal()
                snvLoad += dndsMutationalLoad.snvTotal()
            }
            return CohortLoad(dndsSampleLoad.size, indelLoad, snvLoad)
        }
    }
}
