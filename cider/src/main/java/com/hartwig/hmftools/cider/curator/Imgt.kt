package com.hartwig.hmftools.cider.curator

import com.hartwig.hmftools.cider.IgTcrFunctionality
import com.hartwig.hmftools.cider.IgTcrRegion
import htsjdk.samtools.reference.FastaSequenceFile
import java.io.File

data class ImgtGeneAllele(
    val geneName: String,
    val allele: String,
    val species: String,
    val functionality: IgTcrFunctionality,
    val region: IgTcrRegion?,
    val sequenceWithGaps: String,
    val partial: Boolean)
{
    val sequenceWithoutGaps: String get() { return sequenceWithGaps.replace(".", "") }

    val geneAllele: String get() { return "$geneName*$allele" }
}

// https://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html
private fun functionalityFromImgtCode(code: String): IgTcrFunctionality
{
    return when (code)
    {
        "F" -> IgTcrFunctionality.FUNCTIONAL
        "ORF" -> IgTcrFunctionality.ORF
        "P" -> IgTcrFunctionality.PSEUDOGENE
        else -> throw IllegalArgumentException("unrecognised IMGT functionality: $code")
    }
}

// IGH has three (alpha, delta and gamma) or four (epsilon and mu) constant domains (CH1 to CH4)
private fun igTcrRegionFromImgtCode(region: String): IgTcrRegion?
{
    return when (region)
    {
        "V-REGION" -> IgTcrRegion.V_REGION
        "D-REGION" -> IgTcrRegion.D_REGION
        "J-REGION" -> IgTcrRegion.J_REGION
        "C-REGION" -> IgTcrRegion.CONSTANT // IGK / IGL
        "CH1" -> IgTcrRegion.CONSTANT // IGH
        "EX1" -> IgTcrRegion.CONSTANT // TCR
        else -> null
    }
}

/*
    The FASTA header of IMGT/GENE-DB reference sequences is standardized. It contains 15 fields separated by '|':

    1. IMGT/LIGM-DB accession number(s)
    2. IMGT gene and allele name
    3. species
    4. IMGT allele functionality
    5. exon(s), region name(s), or extracted label(s)
    6. start and end positions in the IMGT/LIGM-DB accession number(s)
    7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
    8. codon start, or 'NR' (not relevant) for non coding labels
    9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
    10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
    11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
    12. number of amino acids (AA): this field indicates that the sequence is in amino acids
    13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
    14. partial (if it is)
    15. reverse complementary (if it is)
*/
// >IMGT000128|IGHA1*06|Homo sapiens|F|M|g,1187575..1187786|213 nt|1|+1| | | |213+0=213| | |
// https://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html
private fun parseGeneData(seqName: String, sequenceWithGaps: String): ImgtGeneAllele
{
    val tokens = seqName.split('|')
    require(tokens.size > 10)

    val geneAllele = tokens[1].split('*')
    require(geneAllele.size == 2)

    // we don't distinguish between F, (F) and [F]
    val functionality = functionalityFromImgtCode(tokens[3].replace(Regex("[()\\[\\]]"), ""))

    val partial = tokens[13].contains("partial")

    return ImgtGeneAllele(
        geneName = geneAllele[0], allele = geneAllele[1], species = tokens[2], functionality = functionality,
        region = igTcrRegionFromImgtCode(tokens[4]), sequenceWithGaps = sequenceWithGaps, partial = partial)
}

fun readGeneDataFromFasta(imgtFastaPath: String): List<ImgtGeneAllele>
{
    val alleles: MutableList<ImgtGeneAllele> = ArrayList()
    val imgtFastaFile = FastaSequenceFile(File(imgtFastaPath), false)
    while (true)
    {
        val sequence = imgtFastaFile.nextSequence() ?: break
        alleles.add(parseGeneData(sequence.name, sequence.baseString.uppercase()))
    }
    return alleles
}
