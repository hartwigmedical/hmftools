package com.hartwig.hmftools.extensions.cli.options.enums

import com.hartwig.hmftools.extensions.cli.options.HmfOption
import com.hartwig.hmftools.extensions.cli.options.validators.OptionValidator
import com.hartwig.hmftools.extensions.cli.options.validators.primitives.EnumValidator
import org.apache.commons.cli.Option

data class EnumOption<T : Enum<T>>(override val name: String, private val enumClass: Class<T>) : HmfOption {
    override val description = "Expected values (case insensitive): ${enumClass.enumConstants.joinToString(", ") { it.name }}"
    override val option: Option = Option.builder(name).desc(description).hasArg().build()
    override val validators: List<OptionValidator> = listOf(EnumValidator(enumClass))
}
