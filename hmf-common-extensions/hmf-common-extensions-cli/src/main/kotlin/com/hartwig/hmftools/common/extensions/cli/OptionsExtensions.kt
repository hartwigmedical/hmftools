package com.hartwig.hmftools.common.extensions.cli

import org.apache.commons.cli.*
import kotlin.system.exitProcess

fun Options.printHelpAndExit(cmd: String, errorMessage: String): Nothing {
    val formatter = HelpFormatter()
    formatter.printHelp(100, cmd, cmd, this, errorMessage, true)
    exitProcess(1)
}

fun Options.createCommandLine(cmd: String, args: Array<String>): CommandLine {
    val parser = DefaultParser()
    return try {
        parser.parse(this, args)
    } catch (parseException: ParseException) {
        this.printHelpAndExit(cmd, "\n${parseException.message}")
    }
}
