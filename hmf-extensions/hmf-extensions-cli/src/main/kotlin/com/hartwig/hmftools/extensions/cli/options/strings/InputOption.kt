package com.hartwig.hmftools.extensions.cli.options.strings

import com.hartwig.hmftools.extensions.cli.options.HmfOption
import com.hartwig.hmftools.extensions.cli.options.validators.OptionValidator
import org.apache.commons.cli.Option

data class InputOption(override val name: String, override val description: String) : HmfOption {
    override val option: Option = Option.builder(name).desc(description).hasArg().build()
    override val validators: List<OptionValidator> = listOf()
}
