package com.hartwig.hmftools.extensions.cli.options.filesystem

import com.google.common.io.Resources
import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import io.kotlintest.matchers.shouldBe
import io.kotlintest.matchers.shouldThrow
import io.kotlintest.specs.StringSpec
import java.io.IOException

class InputFileOptionTest : StringSpec() {

    private val TEST_FILE_PATH = Resources.getResource("cli/testFile.txt").path
    private val TEST_FILE_OPTION = "test_file"

    init {
        "works for valid optional input file" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(InputFileOption(TEST_FILE_OPTION, "Path to test file."))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf("-$TEST_FILE_OPTION", TEST_FILE_PATH))
            optionsWrapper.validate(cmd)
            cmd.hasOption(TEST_FILE_OPTION) shouldBe true
            cmd.getOptionValue(TEST_FILE_OPTION) shouldBe TEST_FILE_PATH
        }

        "works for missing optional input file" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(InputFileOption(TEST_FILE_OPTION, "Path to test file."))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf())
            optionsWrapper.validate(cmd)
            cmd.hasOption(TEST_FILE_OPTION) shouldBe false
        }

        "throws on missing arg for optional input file" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(InputFileOption(TEST_FILE_OPTION, "Path to test file."))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf("-$TEST_FILE_OPTION", "missing_file.txt"))
            shouldThrow<IOException> { optionsWrapper.validate(cmd) }
        }
    }
}
