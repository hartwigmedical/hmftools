package com.hartwig.hmftools.extensions.cli.options

import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import java.io.IOException

class HmfOptions {
    val options = Options()
    private val optionsMap = mutableMapOf<Option, HmfOption>()

    fun add(hmfOption: HmfOption) {
        if (optionsMap.keys.map { it.opt }.contains(hmfOption.name)) {
            println("Warning: Overriding previously set option ${hmfOption.name}")
        }
        optionsMap[hmfOption.option] = hmfOption
        options.addOption(hmfOption.option)
    }

    fun validate(cmd: CommandLine) {
        val parsedOptions = cmd.options.toSet()
        val errorMessages = parsedOptions.mapNotNull { optionsMap[it] }
                .mapNotNull { option -> option.validators.mapNotNull { it.validate(option, cmd) }.firstOrNull() }
                .map { "   $it" }
        if (errorMessages.isNotEmpty()) {
            throw IOException("Could not validate command line options:\n${errorMessages.joinToString("\n")}\n")
        }
    }
}
