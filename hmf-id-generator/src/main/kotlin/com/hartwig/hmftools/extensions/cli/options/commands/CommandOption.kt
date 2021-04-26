package com.hartwig.hmftools.extensions.cli.options.commands

import com.hartwig.hmftools.extensions.cli.options.HmfOption
import com.hartwig.hmftools.extensions.cli.options.validators.OptionValidator
import com.hartwig.hmftools.extensions.cli.options.validators.primitives.CommandValidator
import org.apache.commons.cli.Option

data class CommandOption(override val name: String, override val description: String) : HmfOption {
    override val option: Option = Option.builder(name).desc(description).hasArg().build()
    override val validators: List<OptionValidator> = listOf(CommandValidator)
}
