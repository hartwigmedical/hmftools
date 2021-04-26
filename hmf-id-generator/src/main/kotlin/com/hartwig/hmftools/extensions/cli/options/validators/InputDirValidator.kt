package com.hartwig.hmftools.extensions.cli.options.validators

import com.hartwig.hmftools.extensions.cli.options.HmfOption
import com.hartwig.hmftools.extensions.cli.options.validators.primitives.FileSystemExistsValidator
import com.hartwig.hmftools.extensions.cli.options.validators.primitives.IsDirectoryValidator
import org.apache.commons.cli.CommandLine

object InputDirValidator : OptionValidator {
    override fun validate(option: HmfOption, cmd: CommandLine): String? {
        return FileSystemExistsValidator.validate(option, cmd) ?: IsDirectoryValidator.validate(option, cmd)
    }
}
