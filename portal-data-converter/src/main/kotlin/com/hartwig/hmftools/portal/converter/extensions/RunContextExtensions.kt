package com.hartwig.hmftools.portal.converter.extensions

import com.hartwig.hmftools.common.context.RunContext
import java.io.File

private const val SOMATIC_VCF_EXTENSION_V3 = "_post_processed_v2.2.vcf.gz"
private const val SOMATIC_VCF_EXTENSION_V4 = "_post_processed.vcf.gz"

fun RunContext.somaticVcfPath(): String? {
    return if (!isSomaticRun) null
    else {
        val setName = "${refSample()}_${tumorSample()}"
        val somaticVariantsFolder = "${runDirectory()}/somaticVariants/$setName"
        val somaticVcfV4 = File("$somaticVariantsFolder/$setName$SOMATIC_VCF_EXTENSION_V4")
        val somaticVcfV3 = File("$somaticVariantsFolder/$setName$SOMATIC_VCF_EXTENSION_V3")
        when {
            somaticVcfV4.exists() -> somaticVcfV4.path
            somaticVcfV3.exists() -> somaticVcfV3.path
            else                  -> null
        }
    }
}
