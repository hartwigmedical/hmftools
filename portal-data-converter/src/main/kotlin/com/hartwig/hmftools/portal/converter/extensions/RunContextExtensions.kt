package com.hartwig.hmftools.portal.converter.extensions

import com.hartwig.hmftools.common.context.RunContext
import java.io.File

fun RunContext.somaticVcfPath(): String? {
    if (!isSomaticRun) return null
    else {
        val somaticVariantsFolder = "${runDirectory()}/somaticVariants/${refSample()}_${tumorSample()}"
        val somaticVcf = File("$somaticVariantsFolder/${refSample()}_${tumorSample()}_post_processed.vcf")
        return if (!somaticVcf.exists()) null
        else return somaticVcf.path
    }
}
