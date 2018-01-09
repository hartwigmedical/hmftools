package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.extensions.cli.createCommandLine
import com.hartwig.hmftools.idgenerator.extensions.writeHmfIds
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

fun main(args: Array<String>) {
    val options = createOptions()
    val cmd = options.createCommandLine("hmf-id-create", args)
    val idGenerator = IdGenerator(cmd.getOptionValue(PASSWORD))
    val patients = Files.lines(Paths.get(cmd.getOptionValue(PATIENT_FILE))).map { it.trim() }.filter { it.isNotEmpty() }.collect(Collectors.toList())
    val ids = idGenerator.generateIds(patients)
    val outputFile = File(cmd.getOptionValue(OUTPUT_FILE))
    outputFile.writeHmfIds(ids)
}

private fun createOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(PASSWORD).required().hasArg().desc("password").build())
    options.addOption(Option.builder(PATIENT_FILE).required().hasArg().desc("file containing list of patients to generate ids for").build())
    options.addOption(Option.builder(OUTPUT_FILE).required().hasArg().desc("output file location").build())
    return options
}
