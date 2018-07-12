package com.hartwig.hmftools.extensions.cli.options.validators.primitives

import com.hartwig.hmftools.extensions.cli.options.HmfOption
import com.hartwig.hmftools.extensions.cli.options.validators.OptionValidator
import org.apache.commons.cli.CommandLine

data class EnumValidator<T : Enum<T>>(private val enumClass: Class<T>) : OptionValidator {
    override fun validate(option: HmfOption, cmd: CommandLine): String? {
        if (!cmd.hasOption(option.name)) return null
        val optionValue = cmd.getOptionValue(option.name).toLowerCase()
        val valueSet = enumClass.enumConstants.map { it.name.toLowerCase() }.toSet()
        if (!valueSet.contains(optionValue)) return "-${option.name}: $optionValue is not a valid choice."
        return null
    }
}
