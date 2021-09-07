package com.hartwig.hmftools.extensions

import com.hartwig.hmftools.bedpe.Breakend
import com.hartwig.hmftools.gripss.ContigComparator
import com.hartwig.hmftools.gripss.VariantType
import htsjdk.variant.variantcontext.VariantContext
import org.apache.logging.log4j.LogManager
import org.apache.logging.log4j.util.Strings
import java.lang.Math.abs
import java.util.*
import kotlin.math.max

const val CIPOS = "CIPOS"
const val CIRPOS = "CIRPOS"
const val BEALN = "BEALN"
const val REPEAT_MASKER_REPEAT_CLASS = "INSRMRC"
const val REPEAT_MASKER_REPEAT_TYPE = "INSRMRT"



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

fun VariantContext.realigned(): Boolean {
    return this.getAttributeAsBoolean("REALIGN", false)
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

fun VariantContext.inexactHomologyStart(): Int {
    if (!this.hasAttribute("IHOMPOS")) {
        return 0
    }

    val (iHom1, iHom2) = this.getAttributeAsIntList("IHOMPOS", 0)
    return abs(iHom1)
}

fun VariantContext.inexactHomologyEnd(): Int {
    if (!this.hasAttribute("IHOMPOS")) {
        return 0
    }

    val (iHom1, iHom2) = this.getAttributeAsIntList("IHOMPOS", 0)
    return abs(iHom2)
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

fun VariantContext.hasViralSequenceAlignment(comparator: ContigComparator): Boolean {
    val potentialAlignments = potentialAlignmentLocations();
    return potentialAlignments.isNotEmpty() && !comparator.isValidContig(potentialAlignments[0].contig.replace("chr", ""))
}

fun VariantContext.potentialAlignmentLocations(): List<Breakend> {
    val logger = LogManager.getLogger(this.javaClass)

    if (this.hasAttribute(BEALN)) {
        val result = mutableListOf<Breakend>()
        val stringBeans = this.getAttributeAsStringList(BEALN, "");
        for (stringBean in stringBeans) {
            try {
                result.add(Breakend.fromBealn(stringBean))
            } catch (e: Exception) {
                logger.warn("Malformed bealn field: $stringBean")
            }
        }

        return result
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
