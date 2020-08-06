package com.hartwig.hmftools.extensions

import com.hartwig.hmftools.bedpe.Breakend
import com.hartwig.hmftools.gripss.ContigComparator
import com.hartwig.hmftools.gripss.VariantType
import htsjdk.variant.variantcontext.VariantContext
import org.apache.logging.log4j.util.Strings
import java.util.*
import kotlin.math.max

const val CIPOS = "CIPOS"
const val CIRPOS = "CIRPOS"

fun VariantContext.breakendAssemblyReadPairs(): Int {
    return this.getAttributeAsInt("BASRP", 0)
}

fun VariantContext.breakendAssemblyReadPairsIsInconsistent(): Boolean {
    // This shouldn't occur but sometimes can because of redefined variants in backport
    return this.getAttributeAsInt("BASRP", 0) == 0
            && this.getAttributeAsInt("BASSR", 0) == 0
            && this.getAttributeAsDouble("BAQ", 0.00) > 0
}

fun VariantContext.strandBias(): Double {
    val strandBias = this.getAttributeAsDouble("SB", 0.5)
    return max(strandBias, 1 - strandBias)
}

fun VariantContext.imprecise(): Boolean {
    return this.getAttributeAsBoolean("IMPRECISE", false)
}

fun VariantContext.homologyLength(): Int {
    return this.getAttributeAsInt("HOMLEN", 0)
}

fun VariantContext.homologySequence(): String {
    return this.getAttributeAsString("HOMSEQ", "");
}

fun VariantContext.inexactHomologyLength(): Int {
    if (!this.hasAttribute("IHOMPOS")) {
        return 0
    }

    val (iHom1, iHom2) = this.getAttributeAsIntList("IHOMPOS", 0)
    return max(0, iHom2 - iHom1)
}

fun VariantContext.assemblies(): List<String> {
    val assemblyCount = this.getAttributeAsInt("AS", 0) + this.getAttributeAsInt("RAS", 0) + this.getAttributeAsInt("CAS", 0)
    if (assemblyCount >= 2 && this.hasAttribute("BEID") && this.hasAttribute("BEIDL")) {
        val beid: List<String> = this.getAttributeAsStringList("BEID", Strings.EMPTY)
        val beidl: List<String> = this.getAttributeAsStringList("BEIDL", Strings.EMPTY)
        return beid.zip(beidl) { x, y -> "$x/$y" }
    }

    return Collections.emptyList()
}

fun VariantContext.toVariantType(): VariantType {
    return VariantType.create(this.contig, this.start, this.alleles[0].displayString, this.alleles[1].displayString)
}

fun VariantContext.mate(): String? {
    if (this.hasAttribute("MATEID")) {
        return this.getAttributeAsString("MATEID", "")
    }

    if (this.hasAttribute("PARID")) {
        return this.getAttributeAsString("PARID", "")
    }

    return null
}

fun VariantContext.hasViralSequenceAlignment(): Boolean {
    val potentialAlignments = potentialAlignmentLocations();
    return potentialAlignments.isNotEmpty() && !ContigComparator.defaultContigs.contains(potentialAlignments[0].contig.replace("chr", ""))
}

fun VariantContext.potentialAlignmentLocations(): List<Breakend> {
    if (this.hasAttribute("BEALN")) {
        return this.getAttributeAsStringList("BEALN", "").map { Breakend.fromBealn(it) }
    }

    return listOf()
}

fun VariantContext.confidenceInterval(): Pair<Int, Int> = this.confidenceInterval(CIPOS)

fun VariantContext.remoteConfidenceInterval(): Pair<Int, Int> = this.confidenceInterval(CIRPOS)

private fun VariantContext.confidenceInterval(attribute: String): Pair<Int, Int> {
    if (this.hasAttribute(attribute)) {
        val (start, end) = this.getAttributeAsIntList(attribute, 0)
        return Pair(start, end)
    }
    return Pair(0, 0)
}
