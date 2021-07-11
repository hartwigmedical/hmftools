package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS
import com.hartwig.hmftools.gripss.GripssConfig.Companion.requiredOption
import com.hartwig.hmftools.gripss.GripssHardFilterApplication.Companion.logger
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFFileReader
import org.apache.commons.cli.*
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.IOException

const val INPUT = "input_vcf"
const val OUTPUT = "output_vcf"

fun main(args: Array<String>) {

    @Throws(ParseException::class)
    fun createCommandLine(args: Array<String>, options: Options): CommandLine {
        val parser: CommandLineParser = DefaultParser()
        return parser.parse(options, args)
    }

    fun createOptions(): Options {
        val options = Options()
        options.addOption(requiredOption(INPUT, "Path to input VCF"))
        options.addOption(requiredOption(OUTPUT, "Path to output VCF"))
        return options
    }

    val options = createOptions()
    try {
        val cmd = createCommandLine(args, options)
        val inputVcf = GripssConfig.requiredFile(cmd, INPUT)
        val outputVcf = cmd.getOptionValue(OUTPUT)
        GripssHardFilterApplication(inputVcf, outputVcf).use { x -> x.run() }
    } catch (e: IOException) {
        logger.warn(e)
    } catch (e: ParseException) {
        logger.warn(e)
        val formatter = HelpFormatter()
        formatter.printHelp("gripss", options)
    }
}

class GripssHardFilterApplication(val input: String, val output: String) : Runnable, AutoCloseable {

    companion object {
        val logger = LogManager.getLogger(this::class.java)
        val allowedFilters = setOf(PASS, PON)
    }

    private val fileReader = VCFFileReader(File(input), false)
    private val fileHeader = fileReader.fileHeader
    private val fileWriter = VariantContextWriterBuilder().setReferenceDictionary(fileHeader.sequenceDictionary).setOutputFile(output).build()
    private val startTime = System.currentTimeMillis()

    override fun run() {
        logger.info("Reading $input")
        logger.info("Writing $output")
        fileWriter.writeHeader(fileHeader)
        for (context in fileReader) {
            if (context.isOkayFilter()) {
                fileWriter.add(context.maintainLinkedBy())
            }
        }
    }

    override fun close() {
        fileReader.close()
        fileWriter.close()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")

    }

    private fun VariantContext.maintainLinkedBy(): VariantContext {
        setDotToEmpty(LOCAL_LINKED_BY)
        setDotToEmpty(REMOTE_LINKED_BY)
        return this
    }

    private fun VariantContext.setDotToEmpty(attribute: String) {
        if (this.commonInfo.hasAttribute(attribute) && this.commonInfo.getAttribute(attribute) == ".") {
            this.commonInfo.putAttribute(attribute, "", true)
        }
    }

    private fun VariantContext.isOkayFilter(): Boolean {
        if (this.isNotFiltered) {
            return true
        }

        return this.filters.let { it.size == 1 && it.iterator().next() in allowedFilters }
    }
}