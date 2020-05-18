package com.hartwig.hmftools.gripss

import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.*


const val MIN_QUAL = "minQual";
const val IMPRECISE = "imprecise";
const val MIN_TUMOR_AF = "minTumorAF";
const val SHORT_SR_NORMAL = "shortSRNormal"
const val LONG_DP_SUPPORT = "longDPSupport"
const val SHORT_SR_SUPPORT = "shortSRSupport"
const val MAX_POLY_G_LENGTH = "maxPolyGLength"
const val SHORT_STRAND_BIAS = "shortStrandBias"
const val MAX_NORMAL_SUPPORT = "maxNormalSupport";
const val MIN_NORMAL_COVERAGE = "minNormalCoverage";

const val MAX_HOM_LENGTH = "maxHomLength"
const val MAX_HOM_LENGTH_SHORT_INV = "maxHomLengthShortInv"
const val MAX_INEXACT_HOM_LENGTH = "maxInexactHomLength"
const val MAX_INEXACT_HOM_LENGTH_SHORT_DEL = "maxInexactHomLengthShortDel"
const val BREAK_END_ASSEMBLY_READ_PAIR = "breakendAssemblyReadPair"

const val PON = "PON"
const val PASS = "PASS"
const val DEDUP = "dedup"
const val MIN_SIZE = "minSize"

const val TAF = "TAF";
const val LOCAL_LINKED_BY = "LOCAL_LINKED_BY";
const val REMOTE_LINKED_BY = "REMOTE_LINKED_BY";

class GripssVCF(outputVCF: String) : AutoCloseable {
    val writer = VariantContextWriterBuilder().setOutputFile(outputVCF).unsetOption(Options.INDEX_ON_THE_FLY).build()


    fun writeHeader(header: VCFHeader) {

        //TODO:
        //        ##FILTER=<ID=small.replacement.fp,Description="Deletion with insertion of the same length that is not a simple inversion.">
        //        ##FILTER=<ID=NO_ASRP,Description="Breakend supported by 0 assembled read pairs">
        //        ##FILTER=<ID=cohortMinSize,Description="Variant is smaller than the minimum event size considered for this cohort">

        header.addMetaDataLine(VCFFilterHeaderLine(DEDUP, "Event is duplicate of another"))
        header.addMetaDataLine(VCFFilterHeaderLine(MIN_SIZE, "Event is too short"))
        header.addMetaDataLine(VCFFilterHeaderLine(BREAK_END_ASSEMBLY_READ_PAIR, "Breakend supported by 0 assembled read pairs"))
        header.addMetaDataLine(VCFFilterHeaderLine(MAX_HOM_LENGTH, "Breakpoint homology length too long"))
        header.addMetaDataLine(VCFFilterHeaderLine(MAX_INEXACT_HOM_LENGTH, "Inexact breakpoint homology length too long"))
        header.addMetaDataLine(VCFFilterHeaderLine(MAX_INEXACT_HOM_LENGTH_SHORT_DEL, "Short deletion that appears to be a ligation artifact"))
        header.addMetaDataLine(VCFFilterHeaderLine(MAX_HOM_LENGTH_SHORT_INV, "Short inversion with significant sequence homology"))
        header.addMetaDataLine(VCFFilterHeaderLine(SHORT_SR_SUPPORT, "Short event not supported by any split reads either directly or via assembly"))
        header.addMetaDataLine(VCFFilterHeaderLine(SHORT_SR_NORMAL, "Short event with split reads support in the normal sample"))
        header.addMetaDataLine(VCFFilterHeaderLine(LONG_DP_SUPPORT, "Large event not supported by any read pairs either directly or via assembly"))
        header.addMetaDataLine(VCFFilterHeaderLine(PON, "Found in panel of normals"))
        header.addMetaDataLine(VCFFilterHeaderLine(MIN_TUMOR_AF, "Variant allele fraction too low"))
        header.addMetaDataLine(VCFFilterHeaderLine(MAX_NORMAL_SUPPORT, "Too many support reads from the normal sample"))
        header.addMetaDataLine(VCFFilterHeaderLine(MIN_NORMAL_COVERAGE, "Insufficient normal coverage to determine somatic status"))
        header.addMetaDataLine(VCFFilterHeaderLine(SHORT_STRAND_BIAS, "Short event with excessive strand bias in split reads/soft clipped reads overlapping breakpoint"))
        header.addMetaDataLine(VCFFilterHeaderLine(MIN_QUAL, "Insufficient quality"))
        header.addMetaDataLine(VCFFilterHeaderLine(MAX_POLY_G_LENGTH, "Single breakend containing long polyC or polyG run. Likely to be an artefact"))
        header.addMetaDataLine(VCFFilterHeaderLine(IMPRECISE, "Imprecise variant"))
        header.addMetaDataLine(VCFFilterHeaderLine(PASS, "Variant passes all filters"))
        header.addMetaDataLine(VCFInfoHeaderLine(TAF, 1, VCFHeaderLineType.Float, "Description"))
        header.addMetaDataLine(VCFInfoHeaderLine(LOCAL_LINKED_BY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Breakend linking information"))
        header.addMetaDataLine(VCFInfoHeaderLine(REMOTE_LINKED_BY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Partner breakend linking information"))

        writer.writeHeader(header)
    }

    fun writeVariant(context: VariantContext) {
        writer.add(context)
    }

    override fun close() {
        println("Closing file")
        writer.close()
    }

}