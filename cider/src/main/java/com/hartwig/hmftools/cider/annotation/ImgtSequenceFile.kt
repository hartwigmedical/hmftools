package com.hartwig.hmftools.cider.annotation

import com.hartwig.hmftools.cider.CiderUtils.getResourceAsFile
import com.hartwig.hmftools.cider.CiderUtils.getResourceAsStream
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import htsjdk.samtools.SAMSequenceDictionary
import htsjdk.samtools.reference.FastaSequenceIndex
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import kotlin.io.path.Path

// Curated FASTA file which contains V/D/J gene alleles.
// Includes the exact sequence from IMGT and maybe some surrounding ref context.
class ImgtSequenceFile(genomeVersion: RefGenomeVersion)
{
    val fastaDict: SAMSequenceDictionary
    val bwamemImgPath: String
    val sequencesByContig: Map<String, Sequence>

    init
    {
        val fastaName = "igtcr_gene.${genomeVersion.identifier()}.fasta"

        fastaDict = ReferenceSequenceFileFactory.loadDictionary(getResourceAsStream("$fastaName.dict"))
        bwamemImgPath = getResourceAsFile("$fastaName.img")

        val fastaIndex = FastaSequenceIndex(getResourceAsStream("$fastaName.fai"))
        // Needs to end in .fasta or IndexedFastaSequenceFile will complain.
        val fastaPath = getResourceAsFile(fastaName, ".fasta")
        val fasta = IndexedFastaSequenceFile(Path(fastaPath), fastaIndex)
        sequencesByContig = fastaDict.sequences
            .associateBy(
                { it.sequenceName},
                { Sequence.fromFasta(it.sequenceName, fasta.getSequence(it.sequenceName).baseString) })
    }

    data class Sequence(
        val geneName: String,
        val allele: String,
        val sequenceWithRef: String,
        val refBefore: Int,
        val refAfter: Int
    )
    {
        init
        {
            require(refBefore + refAfter < sequenceWithRef.length)
        }

        val geneAllele: String get() = "$geneName*$allele"
        val imgtRange: IntRange get() = refBefore until (sequenceWithRef.length - refAfter)
        val fastaLabel: String get() = "$geneName|$allele|$refBefore|$refAfter"

        companion object
        {
            fun fromFasta(label: String, sequence: String): Sequence
            {
                val parts = label.split('|')
                val geneName = parts[0]
                val allele = parts[1]
                val refBefore = parts[2].toInt()
                val refAfter = parts[3].toInt()
                return Sequence(geneName, allele, sequence, refBefore, refAfter)
            }
        }
    }
}
