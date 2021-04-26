package com.hartwig.hmftools.extensions.cli.options.validators.primitives

import com.hartwig.hmftools.extensions.cli.options.HmfOption
import com.hartwig.hmftools.extensions.cli.options.validators.OptionValidator
import org.apache.commons.cli.CommandLine
import java.io.File

object FileSystemExistsValidator : OptionValidator {

    override fun validate(option: HmfOption, cmd: CommandLine): String? {
        if (!cmd.hasOption(option.name)) return null
        val file = File(cmd.getOptionValue(option.name))
        if (!file.exists()) return "-${option.name}: $file does not exist."
        return null
    }
}
