
package com.hartwig.hmftools.cider

/*

import com.hartwig.hmftools.common.samtools.CigarUtils
import htsjdk.samtools.SAMRecord
import org.apache.logging.log4j.LogManager
import org.eclipse.collections.api.collection.ImmutableCollection

// read candidate factory is helper class
// to create VJReadCandidate objects, with selected
// preprocessing to make downstream processing simpler
class VJReadCandidateFactory(private val trimBases: Int)
{
    companion object
    {
        private val sLogger = LogManager.getLogger(VJReadCandidateFactory::class.java)

        private fun numTrailingPolyG(seq: String, sliceEnd: Int) : Int
        {
            for (i in 0 until sliceEnd)
            {
                if (seq[sliceEnd - i - 1] != 'G')
                {
                    return i
                }
            }
            return sliceEnd
        }

        private fun numLeadingPolyC(seq: String, sliceStart: Int) : Int
        {
            for (i in sliceStart until seq.length)
            {
                if (seq[i] != 'C')
                {
                    return i
                }
            }
            return seq.length - sliceStart
        }
    }

    // note: this can be called from different threads
    fun tryCreateReadCandiate(
        read: SAMRecord,
        vjAnchorTemplates: ImmutableCollection<VJAnchorTemplate>,
        templateMatchMethod: VJReadCandidate.MatchMethod,
        useRevComp: Boolean,
        readAnchorStart: Int,
        readAnchorEnd: Int,
        templateLocation: GenomeRegionStrand? = null)
            : VJReadCandidate?
    {
        if (vjAnchorTemplates.isEmpty)
        {
            return null
        }

        val vjAnchorTemplate = vjAnchorTemplates.first()

        // find out the imgt gene type. They should be the same type
        val geneType = vjAnchorTemplate.type

        // check to make sure all the same
        if (vjAnchorTemplates.stream().anyMatch { o: VJAnchorTemplate -> o.type !== geneType })
        {
            sLogger.error("multiple gene types found in same match: {}", vjAnchorTemplates)
            throw RuntimeException("multiple gene types found in same match")
        }

        // this step apply the ploy G trimming, base trimming and also trim away
        // part of the sequence away from the anchor
        val readSlice = toReadSlice(read, useRevComp, geneType.vj, readAnchorStart, readAnchorEnd)

        val templateAnchorAA = vjAnchorTemplate.anchorAminoAcidSequence

        // since we don't actually know whether the aligned part is the anchor sequence, we have to use
        // the soft clip that we think make sense
        val leftSoftClip = CigarUtils.leftSoftClip(read)
        val rightSoftClip = CigarUtils.rightSoftClip(read)

        val readMatch = VJReadCandidate(
            readSlice,
            vjAnchorTemplates, geneType,
            vjAnchorTemplate.anchorSequence,
            templateMatchMethod,
            readAnchorStart - readSlice.sliceStart, // adjust by slice start
            readAnchorEnd - readSlice.sliceStart, // adjust by slice start
            leftSoftClip,
            rightSoftClip)

        // set the similarity score
        readMatch.similarityScore = BlosumSimilarityCalc.calcSimilarityScore(
            geneType.vj, readMatch.templateAnchorSequence, readMatch.anchorSequence)

        val geneNames = vjAnchorTemplates.map({ o: VJAnchorTemplate -> o.name }).distinct().toList()
        sLogger.debug(
            "genes: {} read({}) method({}) anchor range({}-{}) template loc({}) "
                    + "anchor AA({}) template AA({}) similarity({})",
            geneNames, read, templateMatchMethod,
            readMatch.anchorOffsetStart, readMatch.anchorOffsetEnd,
            templateLocation,
            readMatch.anchorAA, templateAnchorAA, readMatch.similarityScore)

        // add a check here to make sure we have not made a mistake somewhere
        if (templateMatchMethod == VJReadCandidate.MatchMethod.BLOSUM && readMatch.similarityScore <= 0)
        {
            sLogger.error("blosum match with -ve similarity score: {}", readMatch)
            throw RuntimeException("blosum match with -ve similarity score: ${readMatch}")
        }

        return readMatch
    }

    // apply trim bases and polyG trimming
    private fun toReadSlice(read: SAMRecord, useReverseComplement: Boolean, vj: VJ,
                            readAnchorStart: Int, readAnchorEnd: Int) : ReadSlice
    {
        // work out the slice start and end
        var sliceStart: Int = trimBases
        var sliceEnd: Int = read.readLength - trimBases

        // now we also want to try poly G tail trimming
        // we want to work out there the tail is.
        // the tail is on the right side and poly G if !read.readNegativeStrandFlag
        // the tail is on the left side and poly C otherwise
        if (read.readNegativeStrandFlag)
        {
            val numCs = numLeadingPolyC(read.readString, sliceStart)
            if (numCs >= CiderConstants.MIN_POLY_G_TRIM_COUNT)
            {
                sLogger.info("read: {}, poly G tail found: {}", read, read.readString)
                sliceStart += numCs + CiderConstants.POLY_G_TRIM_EXTRA_BASE_COUNT
            }
        }
        else
        {
            // ends with poly G, but take trim bases into account
            val numGs = numTrailingPolyG(read.readString, sliceEnd)
            if (numGs >= CiderConstants.MIN_POLY_G_TRIM_COUNT)
            {
                sLogger.info("read: {}, poly G tail found: {}", read, read.readString)
                sliceEnd -= numGs + CiderConstants.POLY_G_TRIM_EXTRA_BASE_COUNT
            }
        }

        if (vj == VJ.V)
        {
            // for V, we layout from the anchor start, left to right
            // we are only interested in what comes after anchor start
            sliceStart = Math.max(readAnchorStart, sliceStart)
        }
        else
        {
            // for J, we layout from the anchor last, right to left
            // we are only interested in what comes before anchor end
            sliceEnd = Math.min(readAnchorEnd, sliceEnd)
        }

        return ReadSlice(read, sliceStart, sliceEnd, useReverseComplement)
    }
}

*/
