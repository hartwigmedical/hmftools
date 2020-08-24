package com.hartwig.hmftools.paddle.cohort

import com.hartwig.hmftools.common.drivercatalog.dnds.DndsMutationalLoad
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsMutationalLoadFile

data class TumorMutationalLoad(val indel: Int, val snv: Int) {
    operator fun plus(other: TumorMutationalLoad): TumorMutationalLoad {
        return TumorMutationalLoad(indel + other.indel, snv + other.snv)
    }
}

data class TumorMutationalLoadSample(val sample: String, val biallelicLoad: TumorMutationalLoad, val nonBiallelicLoad: TumorMutationalLoad, val totalLoad: TumorMutationalLoad) {

    companion object {
        fun fromFile(file: String): List<TumorMutationalLoadSample> {
            return DndsMutationalLoadFile.read(file).map { TumorMutationalLoadSample(it) }
        }

        operator fun invoke(load: DndsMutationalLoad): TumorMutationalLoadSample {
            val biallelicLoad = TumorMutationalLoad(load.indelBiallelic(), load.snvBiallelic())
            val nonBiallelicLoad = TumorMutationalLoad(load.indelNonBiallelic(), load.snvNonBiallelic())
            val totalLoad = TumorMutationalLoad(load.indelTotal(), load.snvTotal())

            return TumorMutationalLoadSample(load.sampleId(), biallelicLoad, nonBiallelicLoad, totalLoad)
        }


    }

    operator fun plus(other: TumorMutationalLoadSample): TumorMutationalLoadSample {
        return TumorMutationalLoadSample("Total", biallelicLoad + other.biallelicLoad, nonBiallelicLoad + other.nonBiallelicLoad, totalLoad + other.totalLoad)
    }

}
