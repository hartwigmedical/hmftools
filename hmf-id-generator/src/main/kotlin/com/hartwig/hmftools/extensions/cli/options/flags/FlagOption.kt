package com.hartwig.hmftools.extensions.cli.options.flags

import com.hartwig.hmftools.extensions.cli.options.HmfOption
import com.hartwig.hmftools.extensions.cli.options.validators.OptionValidator
import org.apache.commons.cli.Option

data class FlagOption(override val name: String, override val description: String) : HmfOption {

    override val option: Option = Option.builder(name).desc(description).build()
    override val validators: List<OptionValidator> = listOf()
}
