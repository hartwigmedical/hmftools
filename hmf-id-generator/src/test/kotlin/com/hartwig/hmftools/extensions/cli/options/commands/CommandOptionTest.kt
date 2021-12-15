package com.hartwig.hmftools.extensions.cli.options.commands

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import io.kotest.matchers.shouldBe
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import java.io.IOException

class CommandOptionTest : StringSpec() {

    private val COMMAND = "more"
    private val MISSING_COMMAND = "this_should_be_missing"
    private val CMD_OPTION = "cmd"

    init {
        "works for valid command" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(CommandOption(CMD_OPTION, "Path to test dir."))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf("-$CMD_OPTION", COMMAND))
            optionsWrapper.validate(cmd)
            cmd.hasOption(CMD_OPTION) shouldBe true
            cmd.getOptionValue(CMD_OPTION) shouldBe COMMAND
        }

        "works for missing optional command" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(CommandOption(CMD_OPTION, "Path to test dir."))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf())
            optionsWrapper.validate(cmd)
            cmd.hasOption(CMD_OPTION) shouldBe false
        }

        "throws on missing arg for optional input dir" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(CommandOption(CMD_OPTION, "Path to test dir."))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf("-$CMD_OPTION", MISSING_COMMAND))
            shouldThrow<IOException> { optionsWrapper.validate(cmd) }
        }
    }
}
