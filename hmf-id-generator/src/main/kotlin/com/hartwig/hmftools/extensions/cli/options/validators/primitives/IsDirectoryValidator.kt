package com.hartwig.hmftools.extensions.cli.options.validators.primitives

import com.hartwig.hmftools.extensions.cli.options.HmfOption
import com.hartwig.hmftools.extensions.cli.options.validators.OptionValidator
import org.apache.commons.cli.CommandLine
import java.io.File

object IsDirectoryValidator : OptionValidator {

    override fun validate(option: HmfOption, cmd: CommandLine): String? {
        if (!cmd.hasOption(option.name)) return null
        val dir = File(cmd.getOptionValue(option.name))
        if (!dir.isDirectory) return "-${option.name}: $dir is not a directory."
        return null
    }
}
