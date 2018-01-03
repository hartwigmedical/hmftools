package com.hartwig.hmftools.idgenerator.extensions

import org.apache.commons.cli.*
import org.apache.logging.log4j.LogManager
import kotlin.system.exitProcess

private val logger = LogManager.getLogger("OptionsExtensions")

fun Options.printHelpAndExit(cmd: String): Nothing {
    val formatter = HelpFormatter()
    formatter.printHelp(cmd, cmd, this, "", true)
    exitProcess(1)
}

fun Options.createCommandLine(cmd: String, args: Array<String>): CommandLine {
    val parser = DefaultParser()
    return try {
        parser.parse(this, args)
    } catch (parseException: ParseException) {
        logger.error(parseException.message)
        this.printHelpAndExit(cmd)
    }
}
