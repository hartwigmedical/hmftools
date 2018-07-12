package com.hartwig.hmftools.extensions.cli

import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import org.apache.commons.cli.*
import java.io.IOException
import kotlin.system.exitProcess

fun Options.printHelpAndExit(cmdName: String, errorMessage: String): Nothing {
    val formatter = HelpFormatter()
    formatter.printHelp(200, cmdName, cmdName, this, errorMessage, true)
    exitProcess(1)
}

fun Options.createCommandLine(cmdName: String, args: Array<String>): CommandLine {
    val parser = DefaultParser()
    return try {
        parser.parse(this, args)
    } catch (parseException: ParseException) {
        this.printHelpAndExit(cmdName, "\n$parseException\n")
    }
}

fun Options.createRunModeCommandLine(cmdName: String, args: Array<String>): CommandLine {
    val parser = DefaultParser()
    return try {
        parser.parse(this, args, true)
    } catch (parseException: ParseException) {
        this.printHelpAndExit(cmdName, "\n$parseException\n")
    }
}

fun HmfOptions.createCommandLine(cmdName: String, args: Array<String>): CommandLine {
    val commandLine = this.options.createCommandLine(cmdName, args)
    return try {
        this.validate(commandLine)
        commandLine
    } catch (ioException: IOException) {
        this.options.printHelpAndExit(cmdName, "\n$ioException\n")
    }
}
