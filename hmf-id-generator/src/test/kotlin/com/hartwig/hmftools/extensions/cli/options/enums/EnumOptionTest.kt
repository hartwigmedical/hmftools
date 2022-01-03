package com.hartwig.hmftools.extensions.cli.options.enums

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import io.kotest.matchers.shouldBe
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.StringSpec
import java.io.IOException

class EnumOptionTest : StringSpec() {

    private val MODE = "mode"

    enum class Modes {
        Mode1, @Suppress("unused")
        Mode2
    }

    init {
        "works for valid enum option" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(EnumOption(MODE, Modes::class.java))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf("-$MODE", Modes.Mode1.name))
            optionsWrapper.validate(cmd)
            cmd.hasOption(MODE) shouldBe true
            cmd.getOptionValue(MODE) shouldBe Modes.Mode1.name
        }

        "works for missing enum option" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(EnumOption(MODE, Modes::class.java))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf())
            optionsWrapper.validate(cmd)
            cmd.hasOption(MODE) shouldBe false
        }

        "throws on wrong arg for enum option" {
            val optionsWrapper = HmfOptions()
            optionsWrapper.add(EnumOption(MODE, Modes::class.java))
            val cmd = optionsWrapper.options.createCommandLine("test", arrayOf("-$MODE", "mode3"))
            shouldThrow<IOException> { optionsWrapper.validate(cmd) }
        }
    }
}
