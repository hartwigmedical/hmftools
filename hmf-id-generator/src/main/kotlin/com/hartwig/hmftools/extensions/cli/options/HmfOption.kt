package com.hartwig.hmftools.extensions.cli.options

import com.hartwig.hmftools.extensions.cli.options.validators.OptionValidator
import org.apache.commons.cli.Option

interface HmfOption {

    val name: String
    val description: String
    val option: Option
    val validators: List<OptionValidator>
}
