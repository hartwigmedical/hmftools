package com.hartwig.hmftools.extensions.cli.options.validators.primitives

import com.hartwig.hmftools.extensions.cli.options.HmfOption
import com.hartwig.hmftools.extensions.cli.options.validators.OptionValidator
import org.apache.commons.cli.CommandLine
import org.openqa.selenium.os.ExecutableFinder

object CommandValidator : OptionValidator {
    override fun validate(option: HmfOption, cmd: CommandLine): String? {
        if (!cmd.hasOption(option.name)) return null
        val command = cmd.getOptionValue(option.name)
        if (ExecutableFinder().find(command) == null) return "Command $command (passed as param for -${option.name} option) could not be found."
        return null
    }
}
