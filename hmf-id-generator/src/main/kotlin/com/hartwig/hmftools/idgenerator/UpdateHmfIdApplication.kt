package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.idgenerator.extensions.readOldIds
import com.hartwig.hmftools.idgenerator.extensions.writeHmfIds
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import java.io.File
import java.nio.file.Files
import java.nio.file.Paths
import java.util.stream.Collectors

fun main(args: Array<String>) {
    val options = updateOptions()
    val cmd = options.createCommandLine("hmf-id-update", args)
    val patientFile = cmd.getOptionValue(PATIENT_FILE)
    val generator = IdGenerator(cmd.getOptionValue(PASSWORD))
    val hmfIdsFile = File(cmd.getOptionValue(HMF_IDS_FILE))
    val patients = Files.lines(Paths.get(patientFile)).map { it.trim() }.filter { it.isNotEmpty() }.collect(Collectors.toList())
    val oldIds = hmfIdsFile.readOldIds()
    val newIds = generator.updateIds(cmd.getOptionValue(OLD_PASSWORD), patients, oldIds)
    val outputFile = File(cmd.getOptionValue(OUTPUT_FILE))
    outputFile.writeHmfIds(newIds)
}

private fun updateOptions(): Options {
    val options = Options()
    options.addOption(Option.builder(PASSWORD).required().hasArg().desc("password").build())
    options.addOption(Option.builder(OLD_PASSWORD).required().hasArg().desc("password used to generate hashes in HMF ids file").build())
    options.addOption(Option.builder(PATIENT_FILE).required().hasArg().desc("file containing list of patients to generate ids for").build())
    options.addOption(Option.builder(HMF_IDS_FILE).required().hasArg().desc("file containing current mapping of patients to HMF ids").build())
    options.addOption(Option.builder(OUTPUT_FILE).required().hasArg().desc("output file location").build())
    return options
}
