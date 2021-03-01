package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.common.gripss.GripssFilters.*
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.*


const val PON = "PON"
const val IMPRECISE = "imprecise"
const val SHORT_SR_NORMAL = "shortSRNormalSupport"
const val SHORT_SR_SUPPORT = "shortSRTumorSupport"
const val SHORT_STRAND_BIAS = "shortStrandBias"
const val SHORT_DEL_INS_ARTIFACT = "smallDelInsertionArtifact"

const val MAX_POLY_G_LENGTH = "maxPolyGLength"
const val MAX_POLY_A_HOM_LENGTH = "maxPolyAHomLength"

const val MAX_NORMAL_RELATIVE_SUPPORT = "maxNormalRelativeSupport"
const val MIN_NORMAL_COVERAGE = "minNormalCoverage"
const val DISCORDANT_PAIR_SUPPORT = "discordantPairSupport"

const val MAX_HOM_LENGTH_SHORT_INV = "maxHomLengthShortInv"
const val MAX_INEXACT_HOM_LENGTH_SHORT_DEL = "maxInexactHomLengthShortDel"

const val PASS = "PASS"
const val MIN_LENGTH = "minLength"

const val EVENTTYPE = "EVENTTYPE"
const val TAF = "TAF"
const val ALT_PATH = "ALTP"
const val HOTSPOT = "HOTSPOT"
const val REALIGN = "REALIGN"
const val LOCAL_LINKED_BY = "LOCAL_LINKED_BY"
const val REMOTE_LINKED_BY = "REMOTE_LINKED_BY"

class GripssVCF(outputVCF: String, dictionary: SAMSequenceDictionary) : AutoCloseable {
    val writer = VariantContextWriterBuilder().setReferenceDictionary(dictionary).setOutputFile(outputVCF).build()

    fun writeHeader(version: String, header: VCFHeader, samples: List<String>) {

        val metaData = mutableSetOf<VCFHeaderLine>()
        metaData.addAll(header.metaDataInInputOrder)

        metaData.add(VCFHeaderLine("gripssVersion", version))
        metaData.add(VCFFilterHeaderLine(DEDUP, "Event is duplicate of another"))
        metaData.add(VCFFilterHeaderLine(MIN_LENGTH, "Event is too short"))
        metaData.add(VCFFilterHeaderLine(MAX_INEXACT_HOM_LENGTH_SHORT_DEL, "Short deletion that appears to be a ligation artifact"))
        metaData.add(VCFFilterHeaderLine(MAX_HOM_LENGTH_SHORT_INV, "Short inversion with significant sequence homology"))
        metaData.add(VCFFilterHeaderLine(SHORT_SR_SUPPORT, "Short event not supported by any split reads either directly or via assembly"))
        metaData.add(VCFFilterHeaderLine(SHORT_SR_NORMAL, "Short event with split reads support in the normal sample"))
        metaData.add(VCFFilterHeaderLine(DISCORDANT_PAIR_SUPPORT, "Large event not supported by any read pairs either directly or via assembly"))
        metaData.add(VCFFilterHeaderLine(PON, "Found in panel of normals"))
        metaData.add(VCFFilterHeaderLine(MIN_TUMOR_AF, "Variant allele fraction too low"))
        metaData.add(VCFFilterHeaderLine(MAX_NORMAL_RELATIVE_SUPPORT, "Too many support reads from the normal sample relative to the tumor"))
        metaData.add(VCFFilterHeaderLine(MIN_NORMAL_COVERAGE, "Insufficient normal coverage to determine somatic status"))
        metaData.add(VCFFilterHeaderLine(SHORT_STRAND_BIAS, "Short event with excessive strand bias in split reads/soft clipped reads overlapping breakpoint"))
        metaData.add(VCFFilterHeaderLine(MIN_QUAL, "Insufficient quality"))
        metaData.add(VCFFilterHeaderLine(MAX_POLY_G_LENGTH, "Single breakend containing long polyC or polyG run. Likely to be an artifact"))
        metaData.add(VCFFilterHeaderLine(MAX_POLY_A_HOM_LENGTH, "Homology containing long polyA or polyT run"))
        metaData.add(VCFFilterHeaderLine(IMPRECISE, "Imprecise variant"))
        metaData.add(VCFFilterHeaderLine(PASS, "Variant passes all filters"))
        metaData.add(VCFFilterHeaderLine(SHORT_DEL_INS_ARTIFACT, "Filter any short DEL where the insert sequence length + 1 = deletion length. This is a known GRIDSS artifact."))

        metaData.add(VCFInfoHeaderLine(EVENTTYPE, 1, VCFHeaderLineType.String, "Event type [${EventType.values().joinToString { it.toString()} }]"))
        metaData.add(VCFInfoHeaderLine(TAF, 1, VCFHeaderLineType.Float, "Tumor allelic frequency (fragment support / total support)"))
        metaData.add(VCFInfoHeaderLine(ALT_PATH, 1, VCFHeaderLineType.String, "Alternate path"))
        metaData.add(VCFInfoHeaderLine(LOCAL_LINKED_BY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Breakend linking information"))
        metaData.add(VCFInfoHeaderLine(REMOTE_LINKED_BY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Partner breakend linking information"))
        metaData.add(VCFInfoHeaderLine(HOTSPOT, 1, VCFHeaderLineType.Flag, "Variant is a hotspot"))
        metaData.add(VCFInfoHeaderLine(REALIGN, 1, VCFHeaderLineType.Flag, "Variant was realigned"))

        writer.writeHeader(VCFHeader(metaData, samples))
    }

    fun writeVariant(context: VariantContext) {
        writer.add(context)
    }

    override fun close() {
        writer.close()
    }
}