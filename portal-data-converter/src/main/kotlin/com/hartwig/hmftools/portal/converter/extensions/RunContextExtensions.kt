package com.hartwig.hmftools.portal.converter.extensions

import com.hartwig.hmftools.common.context.RunContext
import java.io.File

private const val SOMATIC_VCF_EXTENSION = "_post_processed.vcf.gz"

fun RunContext.somaticVcfPath(): String? {
    return if (!isSomaticRun) null
    else {
        val setName = "${refSample()}_${tumorSample()}"
        val somaticVariantsFolder = "${runDirectory()}/somaticVariants/$setName"
        val somaticVcf = File("$somaticVariantsFolder/$setName$SOMATIC_VCF_EXTENSION")
        when {
            somaticVcf.exists() -> somaticVcf.path
            else                  -> null
        }
    }
}
