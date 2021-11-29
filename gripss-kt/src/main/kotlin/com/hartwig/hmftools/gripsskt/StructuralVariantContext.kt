package com.hartwig.hmftools.gripsskt

import com.hartwig.hmftools.gripsskt.store.Breakend
import com.hartwig.hmftools.gripsskt.store.Breakpoint
import com.hartwig.hmftools.common.gripss.GripssFilters.MIN_QUAL
import com.hartwig.hmftools.common.gripss.GripssFilters.MIN_TUMOR_AF
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.samtools.reference.ReferenceSequence
import htsjdk.samtools.util.Interval
import htsjdk.samtools.util.Locatable
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder

const val SHORT_CALLING_SIZE = 1000
const val SHORT_RESCUE_SIZE = 10000

typealias Cipos = Pair<Int, Int>

class StructuralVariantContext(val context: VariantContext, private val normalOrdinal: Int = 0, private val tumorOrdinal: Int = 1) {
    companion object {
        private val polyGInsert = "G".repeat(16)
        private val polyCInsert = "C".repeat(16)

        private val polyAHomology = "A".repeat(7)
        private val polyTHomology = "T".repeat(7)

        internal val ReferenceSequence.unambiguousNucleotides: String
            get() = this.baseString.unambiguousNucleotides

        internal val String.unambiguousNucleotides: String
            get() = this.splitToSequence('R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', ignoreCase = true, limit = 0).joinToString(separator = "N")

    }

    private val remoteConfidenceInterval = context.remoteConfidenceInterval()
    private val normalGenotype = if (normalOrdinal != -1) context.getGenotype(normalOrdinal) else null
    private val tumorGenotype = context.getGenotype(tumorOrdinal)

    val contig = context.contig!!
    val imprecise = context.imprecise()
    val precise = !imprecise

    val variantType = context.toVariantType()
    val isMobileElementInsertion = variantType.isMobileElementInsertion()
    val orientation = variantType.startOrientation
    val isSingle = variantType is Single
    val isTooShortToRescue = when (variantType) {
        is Duplication -> variantType.length < SHORT_RESCUE_SIZE
        is Deletion -> variantType.length < SHORT_RESCUE_SIZE
        is Insertion -> variantType.length < SHORT_RESCUE_SIZE
        else -> false
    }

    val vcfId: String = context.id
    val mateId: String? = context.mate()
    val confidenceInterval = context.confidenceInterval()
    val start = context.start

    val inLowQualRegion = inLowQualRegion()

    val tumorQual = tumorGenotype.qual(isSingle, isMobileElementInsertion)

    val startBreakend: Breakend = Breakend(contig, start + confidenceInterval.first, start + confidenceInterval.second, orientation)
    val endBreakend: Breakend? = (variantType as? Paired)?.let { Breakend(it.otherChromosome, it.otherPosition + remoteConfidenceInterval.first, it.otherPosition + remoteConfidenceInterval.second, it.endOrientation) }
    val breakpoint: Breakpoint? = endBreakend?.let { Breakpoint(startBreakend, it) }
    val minStart = startBreakend.start

    val maxStart = startBreakend.end
    val insertSequenceLength = variantType.insertSequence.length
    val duplicationLength = (variantType as? Duplication)?.let { it.length + 1 } ?: 0

    internal val isShortDup = variantType is Duplication && variantType.length < SHORT_CALLING_SIZE
    internal val isShortDel = variantType is Deletion && variantType.length < SHORT_CALLING_SIZE
    internal val isShortIns = variantType is Insertion && variantType.length < SHORT_CALLING_SIZE
    internal val isShort = isShortDup || isShortDel || isShortIns

    private val tumorAF = tumorGenotype.allelicFrequency(isSingle, isShort)

    fun realignRemote(other: StructuralVariantContext): StructuralVariantContext {
        return realignRemote(other.start, other.confidenceInterval)
    }

    private fun realignRemote(otherPosition: Int, otherCipos: Pair<Int, Int>): StructuralVariantContext {
        val refAllele = this.context.reference
        val mate = variantType as Paired
        val alleles = listOf(refAllele, Allele.create(mate.altString(otherPosition, refAllele.displayString)))

        val variantContextBuilder = VariantContextBuilder(context)
                .alleles(alleles)
                .attribute("CIRPOS", listOf(otherCipos.first, otherCipos.second))
                .make()

        return StructuralVariantContext(variantContextBuilder, normalOrdinal, tumorOrdinal)
    }

    fun realign(refGenome: IndexedFastaSequenceFile, comparator: ContigComparator): StructuralVariantContext {
        if (insertSequenceLength > 0) {
            return this
        }

        if (precise && !isSingle) {
            val (centeredCipos, centeredRemoteCipos) = centreAlignConfidenceIntervals(comparator)
            if (centeredCipos != confidenceInterval || centeredRemoteCipos != remoteConfidenceInterval) {
                return realignPaired(refGenome, centeredCipos)
            }
        }

        if (imprecise && (confidenceInterval.first != 0 || confidenceInterval.second != 0)) {
            return if (isSingle) {
                sideAlignSingle(refGenome)
            } else {
                sideAlignPaired(refGenome)
            }
        }

        return this
    }

    fun centreAlignConfidenceIntervals(comparator: ContigComparator): Pair<Cipos, Cipos> {
        assert(!isSingle)
        val mate = variantType as Paired
        val invertStart = orientation == mate.endOrientation && comparator.compare(contig, start, mate.otherChromosome, mate.otherPosition) > 0
        val invertEnd = orientation == mate.endOrientation && !invertStart
        val centeredCipos = centreAlignConfidenceInterval(invertStart, confidenceInterval)
        val centeredRemoteCipos = centreAlignConfidenceInterval(invertEnd, remoteConfidenceInterval)
        return Pair(centeredCipos, centeredRemoteCipos)
    }

    private fun sideAlignPaired(refGenome: IndexedFastaSequenceFile): StructuralVariantContext {
        val newCipos = sideAlignConfidenceInterval(orientation, confidenceInterval)
        return realignPaired(refGenome, newCipos)
    }

    private fun realignPaired(refGenome: IndexedFastaSequenceFile, newCipos: Cipos): StructuralVariantContext {
        if (newCipos == confidenceInterval) {
            return this
        }

        val newStart = updatedPosition(start, confidenceInterval, newCipos)
        val newRef = refGenome.getSubsequenceAt(contig, newStart.toLong(), newStart.toLong()).unambiguousNucleotides

        val mate = variantType as Paired
        val alleles = listOf(Allele.create(newRef, true), Allele.create(mate.altString(mate.otherPosition, newRef)))

        val variantContextBuilder = VariantContextBuilder(context)
                .start(newStart.toLong())
                .stop(newStart.toLong())
                .alleles(alleles)
                .attribute(REALIGN, true)
                .attribute("CIPOS", listOf(newCipos.first, newCipos.second))
                .make()

        return StructuralVariantContext(variantContextBuilder, normalOrdinal, tumorOrdinal)
    }

    private fun sideAlignSingle(refGenome: IndexedFastaSequenceFile): StructuralVariantContext {
        val newCipos = sideAlignConfidenceInterval(orientation, confidenceInterval)
        val newStart = updatedPosition(start, confidenceInterval, newCipos)
        val newRef = refGenome.getSubsequenceAt(contig, newStart.toLong(), newStart.toLong()).unambiguousNucleotides

        val mate = variantType as Single
        val alleles = listOf(Allele.create(newRef, true), Allele.create(mate.altString(newRef)))

        val variantContextBuilder = VariantContextBuilder(context)
                .start(newStart.toLong())
                .stop(newStart.toLong())
                .alleles(alleles)
                .attribute(REALIGN, true)
                .attribute("CIPOS", listOf(newCipos.first, newCipos.second))
                .make()

        return StructuralVariantContext(variantContextBuilder, normalOrdinal, tumorOrdinal)
    }

    private fun updatedPosition(position: Int, oldCipos: Cipos, newCipos: Cipos): Int {
        return position + oldCipos.first - newCipos.first
    }


    private fun centreAlignConfidenceInterval(invert: Boolean, cipos: Cipos): Cipos {
        val totalRange = cipos.second - cipos.first
        val newCiposStart = -totalRange / 2
        val newCiposEnd = totalRange + newCiposStart
        return if (invert) {
            Pair(-newCiposEnd, -newCiposStart)
        } else {
            Pair(newCiposStart, newCiposEnd)
        }
    }

    private fun sideAlignConfidenceInterval(orientation: Byte, cipos: Cipos): Cipos {
        return if (orientation == 1.toByte()) {
            Pair(0, cipos.second - cipos.first)
        } else {
            Pair(cipos.first - cipos.second, 0)
        }
    }

    fun context(localLink: String, remoteLink: String, altPath: String?, isHotspot: Boolean, filters: Set<String>): VariantContext {
        val genotypesToWrite = mutableListOf(tumorGenotype)
        normalGenotype?.let { x -> genotypesToWrite.add(x) }

        val builder = VariantContextBuilder(context).genotypes(genotypesToWrite).filters()
        builder.log10PError(tumorQual / -10.0)
                .attribute(TAF, tumorAF)
                .attribute(LOCAL_LINKED_BY, localLink)
                .attribute(REMOTE_LINKED_BY, remoteLink)
                .attribute(HOTSPOT, isHotspot)
                .attribute(EVENTTYPE,variantType.eventType)

        altPath?.let { x -> builder.attribute(ALT_PATH, x) }
        filters.forEach { x -> builder.filter(x) }
        if (filters.isEmpty()) {
            builder.filter(PASS)
        }

        return builder.make()
    }

    fun assemblies(): List<String> = context.assemblies()

    fun isHardFilter(config: GripssFilterConfig, isHotspot: Boolean): Boolean {
        if (tumorQualFilter(config.hardMinTumorQual)) {
            return true
        }

        if (isHotspot) {
            return false
        }

        if (normalSupportRelativeFilter(config.hardMaxNormalRelativeSupport)) {
            return true
        }

        return (normalSupportAbsoluteFilter(config.hardMaxNormalAbsoluteSupport) && normalSupportRelativeFilter(config.softMaxNormalRelativeSupport))
    }

    fun softFilters(config: GripssFilterConfig): Set<String> {
        val result = mutableSetOf<String>()

        if (normalCoverageFilter(config.minNormalCoverage)) {
            result.add(MIN_NORMAL_COVERAGE)
        }

        if (normalSupportRelativeFilter(config.softMaxNormalRelativeSupport)) {
            result.add(MAX_NORMAL_RELATIVE_SUPPORT)
        }

        if (allelicFrequencyFilter(config.minTumorAF)) {
            result.add(MIN_TUMOR_AF)
        }

        if (shortDelInsertArtifact()) {
            result.add(SHORT_DEL_INS_ARTIFACT)
        }

        if (qualFilter(config.minQualBreakEnd, config.minQualBreakPoint)) {
            result.add(MIN_QUAL)
        }

        if (impreciseFilter()) {
            result.add(IMPRECISE)
        }

        if (polyGCInsertFilter(config.polyGCRegion)) {
            result.add(MAX_POLY_G_LENGTH)
        }

        if (polyATHomologyFilter()) {
            result.add(MAX_POLY_A_HOM_LENGTH)
        }

        if (homologyLengthFilterShortInversion(config.maxHomLengthShortInversion)) {
            result.add(MAX_HOM_LENGTH_SHORT_INV)
        }

        if (inexactHomologyLengthShortDelFilter(config.maxInexactHomLengthShortDel)) {
            result.add(MAX_INEXACT_HOM_LENGTH_SHORT_DEL)
        }

        if (shortSplitReadTumorFilter()) {
            result.add(SHORT_SR_SUPPORT)
        }

        if (shortSplitReadNormalFilter()) {
            result.add(SHORT_SR_NORMAL)
        }

        if (strandBiasFilter(config.maxShortStrandBias)) {
            result.add(SHORT_STRAND_BIAS)
        }

        if (discordantPairSupportFilter()) {
            result.add(DISCORDANT_PAIR_SUPPORT)
        }

        if(singleStrandBias())
            result.add(SGL_STRAND_BIAS)

        if(singleInsertSequenceMinLength())
            result.add(SGL_INS_SEQ_MIN_LENGTH)

        if (minLengthFilter(config.minLength)) {
            result.add(MIN_LENGTH)
        }

        return result
    }

    fun qualFilter(minQualBreakEnd: Int, minQualBreakPoint: Int): Boolean {
        var minQual = if (isSingle) minQualBreakEnd else minQualBreakPoint

        if(inLowQualRegion)
            minQual = minQual / 2

        return tumorQual < minQual.toDouble()
    }

    fun polyGCInsertFilter(polyGRegion: Locatable): Boolean {
        return if (isSingle) {
            variantType.insertSequence.contains(polyGInsert) || variantType.insertSequence.contains(polyCInsert) or polyGRegion.contains(context)
        } else {
            polyGRegion.contains(context) || polyGRegion.contains((variantType as Paired).let { Interval(it.otherChromosome, it.otherPosition, it.otherPosition) })
        }
    }

    fun polyATHomologyFilter(): Boolean {
        val homseq = context.homologySequence()
        return homseq.contains(polyAHomology) || homseq.contains(polyTHomology);
    }

    fun inexactHomologyLengthShortDelFilter(maxInexactHomLength: Int, minDelLength: Int = 100, maxDelLength: Int = 800): Boolean {
        return variantType is Deletion && variantType.length >= minDelLength && variantType.length <= maxDelLength && context.inexactHomologyLength() > maxInexactHomLength
    }

    fun inexactHomologyStart(): Int { return context.inexactHomologyStart(); }
    fun inexactHomologyEnd(): Int { return context.inexactHomologyEnd(); }

    fun singleStrandBias(): Boolean {
        if(!isSingle)
            return false

        if(isMobileElementInsertion)
            return false

        return context.strandBias() < 0.05 || context.strandBias() > 0.95
    }

    fun singleInsertSequenceMinLength(): Boolean {
        if(!isSingle)
            return false

        if(isMobileElementInsertion)
            return false

        return variantType.insertSequence.length < 16
    }

    fun minLengthFilter(minSize: Int): Boolean {
        return when (variantType) {
            is Deletion -> variantType.length + variantType.insertSequence.length - 1 < minSize
            is Insertion -> variantType.length + variantType.insertSequence.length + 1 < minSize
            is Duplication -> variantType.length + variantType.insertSequence.length < minSize
            else -> false
        }
    }

    fun homologyLengthFilterShortInversion(maxHomLength: Int, maxInversionLength: Int = 40): Boolean {
        return variantType is Inversion && variantType.length <= maxInversionLength && context.homologyLength() > maxHomLength
    }

    fun shortSplitReadTumorFilter(): Boolean {
        val splitRead = tumorGenotype.splitRead()
        val indelCount = tumorGenotype.indelCount()

        return isShort && (splitRead + indelCount) == 0
    }

    fun shortSplitReadNormalFilter(): Boolean {
        if (normalGenotype == null) {
            return false
        }

        val splitRead = normalGenotype.splitRead()
        val indelCount = normalGenotype.indelCount()

        return isShort && splitRead + indelCount > 0
    }

    fun discordantPairSupportFilter(): Boolean {
        if (normalGenotype == null) {
            return false
        }

        if(isSingle)
            return false

        return !isShort && normalGenotype.readPairs() == 0 && normalGenotype.assemblyReadPairs() == 0
                && tumorGenotype.readPairs() == 0 && tumorGenotype.assemblyReadPairs() == 0
    }

    fun impreciseFilter(): Boolean {
        return context.imprecise()
    }

    fun strandBiasFilter(maxShortStrandBias: Double): Boolean {
        return isShort && context.strandBias() > maxShortStrandBias
    }

    fun allelicFrequencyFilter(minTumorAf: Double): Boolean {
        if(isSingle)
            return tumorAF < 0.015 // turn into config in new GRIPSS
        else
            return tumorAF < minTumorAf
    }

    fun shortDelInsertArtifact(): Boolean {
        if (isShortDel) {
            val deletion = variantType as Deletion
            return deletion.length - 1 == deletion.insertSequence.length
        }
        return false
    }

    fun normalCoverageFilter(minNormalCoverage: Int): Boolean {
        if (normalGenotype == null) {
            return false
        }

        val supportingFragments = normalGenotype.fragmentSupport(isSingle)
        val ref = normalGenotype.refSupportRead()
        val refPair = normalGenotype.refSupportReadPair()

        return supportingFragments + ref + refPair < minNormalCoverage
    }

    fun normalSupportAbsoluteFilter(maxNormalAbsoluteSupport: Int): Boolean {
        if (normalGenotype == null) {
            return false
        }

        val normalSupport = normalGenotype.fragmentSupport(isSingle)
        return normalSupport > maxNormalAbsoluteSupport
    }

    fun tumorQualFilter(minTumorQual: Int): Boolean {
        return tumorQual < minTumorQual.toDouble()
    }

    fun normalSupportRelativeFilter(maxNormalRelativeSupport: Double): Boolean {
        if (normalGenotype == null ) {
            return false
        }

        val normalSupport = normalGenotype.fragmentSupport(isSingle)
        val tumorSupport = tumorGenotype.fragmentSupport(isSingle)

        return normalSupport > maxNormalRelativeSupport * tumorSupport
    }

    override fun toString(): String {
        return "${context.id} ${context.contig}:${context.start} QUAL:${tumorQual} Orientation:${variantType.startOrientation} ${context.alleles[0].displayString}  > ${context.alleles[1].displayString}"
    }

    private fun Cipos.invert(): Cipos {
        return Pair(-this.second, -this.first)
    }

    private fun inLowQualRegion(): Boolean
    {
        // ENSG00000122512	PMS2	7	-1	6012870	6048756
        if(contig.equals("7") || contig.equals("chr7"))
        {
            if(start > 6002870 && start < 6058756)
                return true
        }

        return false
    }



}