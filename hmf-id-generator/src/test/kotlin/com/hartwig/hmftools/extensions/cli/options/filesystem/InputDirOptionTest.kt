package com.hartwig.hmftools.extensions.cli.options.filesystem

import com.google.common.io.Resources
import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import io.kotlintest.matchers.shouldBe
import io.kotlintest.matchers.shouldThrow
import io.kotlintest.specs.StringSpec
import java.io.IOException

class InputDirOptionTest : StringSpec() {

    private val TEST_FILE_PATH = Resources.getResource("testFile.txt").path
    private val TEST_DIR_PATH = Resources.getResource("").path
    private val TEST_DIR_OPTION = "test_dir"

    init {
        "works for valid optional input dir" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(InputDirOption(TEST_DIR_OPTION, "Path to test dir."))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf("-$TEST_DIR_OPTION", TEST_DIR_PATH))
            optionsWrapper.validate(cmd)
            cmd.hasOption(TEST_DIR_OPTION) shouldBe true
            cmd.getOptionValue(TEST_DIR_OPTION) shouldBe TEST_DIR_PATH
        }

        "works for missing optional input dir" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(InputDirOption(TEST_DIR_OPTION, "Path to test dir."))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf())
            optionsWrapper.validate(cmd)
            cmd.hasOption(TEST_DIR_OPTION) shouldBe false
        }

        "throws on missing arg for optional input dir" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(InputDirOption(TEST_DIR_OPTION, "Path to test dir."))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf("-$TEST_DIR_OPTION", "missing_dir"))
            shouldThrow<IOException> { optionsWrapper.validate(cmd) }
        }

        "throws on path to file instead of dir" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(InputDirOption(TEST_DIR_OPTION, "Path to test dir."))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf("-$TEST_DIR_OPTION", TEST_FILE_PATH))
            shouldThrow<IOException> { optionsWrapper.validate(cmd) }
        }
    }
}
