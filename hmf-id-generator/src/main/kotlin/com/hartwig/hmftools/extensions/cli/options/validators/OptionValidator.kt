package com.hartwig.hmftools.extensions.cli.options.validators

import com.hartwig.hmftools.extensions.cli.options.HmfOption
import org.apache.commons.cli.CommandLine

interface OptionValidator {

    fun validate(option: HmfOption, cmd: CommandLine): String?
}
