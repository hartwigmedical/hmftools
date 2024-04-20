package com.hartwig.hmftools.cider

import com.google.common.collect.ImmutableCollection
import com.hartwig.hmftools.cider.VJReadCandidate.MatchMethod
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.cider.genes.IgTcrConstantDiversityRegion
import com.hartwig.hmftools.common.genome.region.GenomeRegion
import com.hartwig.hmftools.common.genome.region.GenomeRegions
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.bam.CigarUtils
import com.hartwig.hmftools.common.bam.SamRecordUtils
import com.hartwig.hmftools.common.utils.IntPair
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.util.CoordMath
import htsjdk.samtools.util.SequenceUtil
import org.apache.logging.log4j.LogManager
import java.util.*
import java.util.concurrent.ConcurrentHashMap

class CiderReadScreener(// collect the reads and sort by types
    private val mCiderGeneDatastore: ICiderGeneDatastore,
    private val mAnchorBlosumSearcher: IAnchorBlosumSearcher,
    private val mMaxReadDistanceFromAnchor: Int,
    private val mMaxFragmentLength: Int
)
{
    // this read key would work for supplementary reads also
    data class ReadRecordKey(
        val readName: String,
        val firstOfPair: Boolean,
        val referenceIndex: Int,
        val alignmentStart: Int
    )

    // want a map to make sure we do not process same record twice
    private val mProcessedReadRecords = Collections.newSetFromMap(ConcurrentHashMap<ReadRecordKey, Boolean>())

    // following store the reads we found, must be thread safe
    private val mVjReadCandidates = Collections.synchronizedList(ArrayList<VJReadCandidate>())
    private val mAllMatchedReads = Collections.synchronizedList(ArrayList<SAMRecord>())

    val vjReadCandidates: List<VJReadCandidate>
        get() = mVjReadCandidates
    val allMatchedReads: List<SAMRecord>
        get() = mAllMatchedReads

    // Note: this function is called from multiple threads
    fun asyncProcessSamRecord(samRecord: SAMRecord)
    {
        // see if we already processed this read. We must check using the referenceIndex. In the case of
        // unmapped read, the read mate position will be the one used. However in the later steps we
        // do not actually want to use this.
        val mapped = if (samRecord.referenceIndex == -1) null else GenomeRegions.create(
            samRecord.referenceName,
            samRecord.alignmentStart, samRecord.alignmentEnd)

        val readRecordKey = ReadRecordKey(samRecord.readName, SamRecordUtils.firstInPair(samRecord),
                                            samRecord.referenceIndex, samRecord.alignmentStart)

        if (!mProcessedReadRecords.add(readRecordKey))
        {
            //sLogger.trace("read record {} already processed", readRecordKey);
            return
        }

        var matchFound = false

        // first we try to match by genome region
        if (mapped != null)
        {
            if (!samRecord.readUnmappedFlag)
            {
                matchFound = matchFound or tryMatchByAlignment(samRecord, mapped)
            }
            else
            {
                // if read is not mapped then the mate is mapped to this region instead, we use the mate mapped
                // region to decide which locus (IGH, TRA etc) we need to search for
                // if we do not do this then we will very likely match to wrong type
                matchFound = matchFound or tryMatchUnmappedReadByBlosum(samRecord)
            }
        }
        else
        {
            // if both read and its mate are unmapped then we just try match by blosum
            matchFound = matchFound or tryMatchByBlosum(samRecord)
        }

        if (matchFound)
        {
            mAllMatchedReads.add(samRecord)
        }
    }

    private fun tryMatchByAlignment(samRecord: SAMRecord, mapped: GenomeRegion): Boolean
    {
        // first step we see if the read overlaps with any anchor location
        for (anchorLocation in mCiderGeneDatastore.getVjAnchorGeneLocations())
        {
            val readCandidate = matchesAnchorLocation(samRecord, mapped, anchorLocation, false)
            if (readCandidate != null)
            {
                return true
            }
        }
        // if none overlaps then we try matching with extrapolation
        for (anchorLocation in mCiderGeneDatastore.getVjAnchorGeneLocations())
        {
            val readCandidate = matchesAnchorLocation(samRecord, mapped, anchorLocation, true)
            if (readCandidate != null)
            {
                return true
            }
        }

        val leftSoftClip = CigarUtils.leftSoftClipLength(samRecord)
        val rightSoftClip = CigarUtils.rightSoftClipLength(samRecord)
        if (leftSoftClip != 0 || rightSoftClip != 0)
        {
            for (igTcrConstantRegion in mCiderGeneDatastore.getIgConstantDiversityRegions())
            {
                // now try to match around location of constant regions
                val readCandidate = tryMatchFromConstantRegion(samRecord, mapped, igTcrConstantRegion)
                if (readCandidate != null)
                {
                    return true
                }
            }
        }
        return false
    }

    // what this function does is to see if the anchor location is within this read
    fun matchesAnchorLocation(
        samRecord: SAMRecord, mapped: GenomeRegion,
        anchorLocation: VJAnchorGenomeLocation, allowExtrapolation: Boolean
    ): VJReadCandidate?
    {
        val readLength = samRecord.readLength

        if (allowExtrapolation)
        {
            if (!isRelevantToAnchorLocation(readLength, mapped, anchorLocation, mMaxReadDistanceFromAnchor))
                return null
        }
        else
        {
            // only if mapped
            if (!isMappedToAnchorLocation(mapped, anchorLocation))
                return null
        }
        val anchorLength: Int = anchorLocation.baseLength()
        if (anchorLength != 30)
        {
            throw RuntimeException("unexpected anchor length")
        }
        val readAnchorRange = extrapolateAnchorReadRange(samRecord, anchorLocation) ?: return null
        var readAnchorStart = readAnchorRange.left
        var readAnchorEnd = readAnchorRange.right

        // anchor does not have to overlap with the read at all
        // we apply extrapolation to get reads that do not overlap with anchor
        if (readAnchorStart < samRecord.readLength + mMaxReadDistanceFromAnchor &&
            readAnchorEnd > -mMaxReadDistanceFromAnchor)
        {
            if (anchorLocation.strand === Strand.REVERSE)
            {
                // say if read length is 2, and read anchor start is 1, read anchor end is 2
                // reverse it we want start = 1, end = 2
                val revStart = samRecord.readLength - readAnchorEnd
                readAnchorEnd = samRecord.readLength - readAnchorStart
                readAnchorStart = revStart
            }

            // want to make sure same gene is not included twice
            val genes: ImmutableCollection<VJAnchorTemplate> = mCiderGeneDatastore.getByGeneLocation(anchorLocation.genomeLocation)

            if (!genes.isEmpty())
            {
                return addVjReadCandidate(samRecord, genes, MatchMethod.ALIGN,
                    anchorLocation.strand === Strand.REVERSE,
                    readAnchorStart, readAnchorEnd, anchorLocation.genomeLocation)
            }
        }
        return null
    }

    // This read is unmapped and mate is mapped. We use the mate mapping to find out
    // which Ig/TCR locus this read is near. And only blosum search for anchor with
    // this locus. This is needed to make sure we find the correct locus
    private fun tryMatchUnmappedReadByBlosum(read: SAMRecord): Boolean
    {
        require(read.readPairedFlag)
        require(read.readUnmappedFlag)
        require(!read.mateUnmappedFlag)

        var relevantAnchorLocation: VJAnchorGenomeLocation? = null
        var relevantConstantRegion: IgTcrConstantDiversityRegion? = null

        // look through anchor locations and find ones that we can use
        for (anchorLocation: VJAnchorGenomeLocation in mCiderGeneDatastore.getVjAnchorGeneLocations())
        {
            if (isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, mMaxFragmentLength))
            {
                relevantAnchorLocation = anchorLocation
                break
            }
        }

        if (relevantAnchorLocation == null)
        {
            // try the constant region
            for (constantRegion: IgTcrConstantDiversityRegion in mCiderGeneDatastore.getIgConstantDiversityRegions())
            {
                if (isUnamppedReadRelevantToConstantRegion(read, constantRegion.genomeLocation, mMaxFragmentLength))
                {
                    relevantConstantRegion = constantRegion
                    break
                }
            }

            if (relevantConstantRegion == null)
            {
                return false
            }
        }

        // get the VJ gene type we are interested in
        // we need to match both the V and J side after an anchor region match
        // i.e. if we match a V gene type, we want to also check for the J gene type, and vice versa

        val vjGeneTypes: List<VJGeneType> =

        if (relevantAnchorLocation != null)
        {
            listOf(relevantAnchorLocation.vjGeneType) + relevantAnchorLocation.vjGeneType.pairedVjGeneTypes()
        }
        else
        {
            relevantConstantRegion!!.getCorrespondingVJ()
        }

        if (vjGeneTypes.isEmpty())
        {
            return false
        }

        //val vjGeneTypes = listOf(relevantVjGeneType) + relevantVjGeneType  relaventVjGeneType.pairedVjGeneTypes()

        // we test both reverse and not reverse, it is simpler this way
        for (reverseRead in arrayOf(true, false))
        {
            var readString = read.readString
            if (reverseRead)
            {
                readString = SequenceUtil.reverseComplement(readString)
            }

            val anchorBlosumMatch: AnchorBlosumMatch? = mAnchorBlosumSearcher.searchForAnchor(
                readString,
                vjGeneTypes,
                IAnchorBlosumSearcher.Mode.DISALLOW_NEG_SIMILARITY,
                0,
                read.readLength
            )
            if (anchorBlosumMatch != null && anchorBlosumMatch.similarityScore > 0)
            {
                val readCandidate = addVjReadCandidate(
                    read,
                    anchorBlosumMatch.templateGenes,
                    MatchMethod.BLOSUM,
                    reverseRead,
                    anchorBlosumMatch.anchorStart,
                    anchorBlosumMatch.anchorEnd
                )
                if (readCandidate != null)
                {
                    if (relevantAnchorLocation != null)
                    {
                        sLogger.trace(
                            "read({}) matched from mate mapped({}:{}) near anchor location({})",
                            read, read.mateReferenceName, read.mateAlignmentStart, relevantAnchorLocation
                        )
                    }
                    else
                    {
                        sLogger.trace(
                            "read({}) matched from mate mapped({}:{}) near constant region({})",
                            read, read.mateReferenceName, read.mateAlignmentStart, relevantConstantRegion
                        )
                    }
                    return true
                }
            }
        }

        return false
    }

    private fun tryMatchByBlosum(samRecord: SAMRecord): Boolean
    {
        // in case someone changes this
        assert(Strand.values().size == 2)

        for (strand in Strand.values())
        {
            var readString = samRecord.readString
            if (strand == Strand.REVERSE) readString = SequenceUtil.reverseComplement(readString)

            val anchorBlosumMatch: AnchorBlosumMatch? =
                mAnchorBlosumSearcher.searchForAnchor(readString, IAnchorBlosumSearcher.Mode.DISALLOW_NEG_SIMILARITY)

            if (anchorBlosumMatch != null && anchorBlosumMatch.similarityScore > 0)
            {
                val readCandidate = addVjReadCandidate(
                    samRecord,
                    anchorBlosumMatch.templateGenes,
                    MatchMethod.BLOSUM,
                    strand == Strand.REVERSE,
                    anchorBlosumMatch.anchorStart,
                    anchorBlosumMatch.anchorEnd
                )
                if (readCandidate != null)
                {
                    return true
                }
            }
        }
        return false
    }

    fun tryMatchFromConstantRegion(
        samRecord: SAMRecord, mapped: GenomeRegion, igTcrConstantDiversityRegion: IgTcrConstantDiversityRegion
    ): VJReadCandidate?
    {
        if (!igTcrConstantDiversityRegion.genomeLocation.inPrimaryAssembly)
            return null

        val readLength = samRecord.readLength
        val (chromosome, posStart, posEnd, strand, _) = igTcrConstantDiversityRegion.genomeLocation

        // see if the anchor location is mapped around here
        if (posStart - readLength < mapped.end() &&
            posEnd + readLength > mapped.start() && chromosome == mapped.chromosome()
        )
        {
            // this is mapped to an IG constant region
            // we want to then see if there is any clipping going on
            // only with soft clip can we say this is actually potentially a J region
            // If these are in positive strand:
            // VVVVV-DDDD-JJ-CCCC
            // If the read is mapped to the constant CCCC region, then the left soft clip might
            // contain J
            // in negative strand, it would be right soft clip
            val anchorBlosumMatch: AnchorBlosumMatch?

            if (strand == Strand.FORWARD)
            {
                val leftSoftClip = CigarUtils.leftSoftClipLength(samRecord)
                if (leftSoftClip == 0) return null

                // now try to find an anchor here
                anchorBlosumMatch = mAnchorBlosumSearcher.searchForAnchor(
                    samRecord.readString,
                    igTcrConstantDiversityRegion.getCorrespondingVJ(),
                    IAnchorBlosumSearcher.Mode.DISALLOW_NEG_SIMILARITY,
                    0,
                    leftSoftClip)
            }
            else
            {
                // reverse strand
                val rightSoftClip = CigarUtils.rightSoftClipLength(samRecord)
                if (rightSoftClip == 0) return null

                val reverseCompSeq = SequenceUtil.reverseComplement(samRecord.readString)

                anchorBlosumMatch =  mAnchorBlosumSearcher.searchForAnchor(
                    reverseCompSeq,
                    igTcrConstantDiversityRegion.getCorrespondingVJ(),
                    IAnchorBlosumSearcher.Mode.DISALLOW_NEG_SIMILARITY,
                    0, rightSoftClip)
            }
            if (anchorBlosumMatch != null && anchorBlosumMatch.similarityScore > 0)
            {
                sLogger.trace("read({}) matched from constant region({}) using blossom", samRecord, igTcrConstantDiversityRegion)

                return addVjReadCandidate(
                    samRecord,
                    anchorBlosumMatch.templateGenes,
                    MatchMethod.BLOSUM,
                    strand == Strand.REVERSE,
                    anchorBlosumMatch.anchorStart,
                    anchorBlosumMatch.anchorEnd
                )
            }
        }
        return null
    }

    private fun addVjReadCandidate(
        samRecord: SAMRecord,
        vjAnchorTemplates: ImmutableCollection<VJAnchorTemplate>,
        templateMatchMethod: MatchMethod,
        useRevComp: Boolean,
        readAnchorStart: Int,
        readAnchorEnd: Int,
        templateLocation: GenomicLocation? = null)
    : VJReadCandidate?
    {
        if (vjAnchorTemplates.isEmpty())
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
        val templateAnchorAA = vjAnchorTemplate.anchorAminoAcidSequence

        // since we don't actually know whether the aligned part is the anchor sequence, we have to use
        // the soft clip that we think make sense
        val leftSoftClip = CigarUtils.leftSoftClipLength(samRecord)
        val rightSoftClip = CigarUtils.rightSoftClipLength(samRecord)

        val readMatch = VJReadCandidate(
            samRecord, vjAnchorTemplates, geneType,
            vjAnchorTemplate.anchorSequence,
            templateMatchMethod, useRevComp,
            readAnchorStart, readAnchorEnd, leftSoftClip, rightSoftClip)

        // set the similarity score
        if (templateMatchMethod == MatchMethod.BLOSUM)
        {
            readMatch.similarityScore = BlosumSimilarityCalc.calcSimilarityScore(
                geneType.vj, readMatch.templateAnchorSequence, readMatch.anchorSequence
            )
        }

        val geneNames = vjAnchorTemplates.map({ o: VJAnchorTemplate -> o.geneName }).distinct().toList()
        sLogger.trace(
            "genes: {} read({}) method({}) anchor range({}-{}) template loc({}) "
                    + "anchor AA({}) template AA({}) similarity({})",
            geneNames, samRecord, templateMatchMethod,
            readMatch.anchorOffsetStart, readMatch.anchorOffsetEnd,
            templateLocation,
            readMatch.anchorAA, templateAnchorAA, readMatch.similarityScore)

        // add a check here to make sure we have not made a mistake somewhere
        if (templateMatchMethod == MatchMethod.BLOSUM && readMatch.similarityScore <= 0)
        {
            sLogger.error("blosum match with -ve similarity score: {}", readMatch)
            throw RuntimeException("blosum match with -ve similarity score: ${readMatch}")
        }

        // add it to list
        mVjReadCandidates.add(readMatch)
        return readMatch
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(CiderReadScreener::class.java)

        fun isMappedToAnchorLocation(mapped: GenomeRegion, anchorLocation: VJAnchorGenomeLocation): Boolean
        {
            if (anchorLocation.chromosome != mapped.chromosome()) return false
            return anchorLocation.start <= mapped.end() && mapped.start() <= anchorLocation.end
        }

        // we are looking for read that can be extrapolated to the anchor location even if they do not overlap
        fun isRelevantToAnchorLocation(readLength: Int, mapped: GenomeRegion,
                                       anchorLocation: VJAnchorGenomeLocation,
                                       maxReadDistanceFromAnchor: Int): Boolean
        {
            if (anchorLocation.chromosome != mapped.chromosome()) return false

            // for V we only allow reads that are mapped upstream
            // for J we only allow reads that are mapped downstream
            //
            // >------------V--------D----------J---------> positive strand
            //             ======           ======
            //
            // or
            //
            // <------------J--------D----------V---------< negative strand
            //             ======           ======
            //
            // reads mapped around the ===== sections are ok
            // we translate it to the genome coord space.
            val anchorLength: Int = anchorLocation.baseLength()
            val isMappedAroundHere: Boolean
            if (anchorLocation.vj == VJ.V && anchorLocation.strand == Strand.FORWARD ||
                anchorLocation.vj == VJ.J && anchorLocation.strand == Strand.REVERSE)
            {
                // here we want the anchor to be downstream
                // we want anchor mapped to higher coord than or equal read
                // we allow anchor to overshoot the mapped region by half the anchor length
                //           |__read__|
                //    |_________________________________________________________|     allowed anchor range
                //     anchor           read length + MaxReadDistanceFromAnchor
                isMappedAroundHere = anchorLocation.start > mapped.start() - anchorLength &&
                        anchorLocation.end < mapped.end() + readLength + maxReadDistanceFromAnchor
            } else
            {
                // here we want the anchor to be upstream
                // we want anchor mapped to lower coord than or equal read
                // we allow anchor to overshoot the mapped region by half the anchor length
                //                                            |__read__|
                // |____________________________________________________________|    allowed anchor range
                //   read length + MaxReadDistanceFromAnchor             anchor
                isMappedAroundHere = anchorLocation.end < mapped.end() + anchorLength &&
                        anchorLocation.start > mapped.start() - readLength - maxReadDistanceFromAnchor
            }
            return isMappedAroundHere
        }

        // 0 based read offset
        // this is similar to SAMRecord getReadPositionAtReferencePosition
        // Two major differences are:
        // 1. this one returns 0 base offset
        // 2. this function will extrapolate into other regions if reference pos is not with any
        //    alignment block. For example, for cigar 20M1000N30M, if we query for a position
        //    10 bases after the 20M block, we will get a read offset of 30 (20M + 10 bases)
        // 3. can return negative or >= read length if the extraplated position is not within read
        fun extrapolateReadOffsetAtRefPosition(record: SAMRecord, referencePos: Int): Int
        {
            // since it is almost impossible to have 3+ alignment blocks speed should not be an issue
            var closesAlignmentBlockDistance = Int.MAX_VALUE
            var readOffset = 0
            for (alignmentBlock in record.alignmentBlocks)
            {
                val blockReferenceEnd = CoordMath.getEnd(alignmentBlock.referenceStart, alignmentBlock.length)
                val blockReferenceStart = alignmentBlock.referenceStart
                val blockDistance: Int

                if (referencePos < blockReferenceStart)
                {
                    blockDistance = blockReferenceStart - referencePos
                }
                else if (referencePos > blockReferenceEnd)
                {
                    blockDistance = referencePos - blockReferenceStart
                }
                else
                {
                    // the anchor reference position is within this block
                    blockDistance = 0
                }
                if (blockDistance < closesAlignmentBlockDistance)
                {
                    closesAlignmentBlockDistance = blockDistance
                    readOffset = alignmentBlock.readStart + referencePos - blockReferenceStart - 1 // make it 0 based
                }
            }

            return readOffset
        }

        fun extrapolateAnchorReadRange(record: SAMRecord, anchorLocation: VJAnchorGenomeLocation): IntPair?
        {
            // we always use the reference position to find it
            val anchorEndReferencePos: Int = anchorLocation.anchorBoundaryReferencePosition()
            val anchorEndReadOffset = extrapolateReadOffsetAtRefPosition(record, anchorEndReferencePos)

            if (anchorEndReferencePos == anchorLocation.start)
            {
                return IntPair(anchorEndReadOffset, anchorEndReadOffset + anchorLocation.baseLength())
            }
            else if (anchorEndReferencePos == anchorLocation.end)
            {
                // make end exclusive
                return IntPair(anchorEndReadOffset - anchorLocation.baseLength() + 1, anchorEndReadOffset + 1)
            }
            else
            {
                throw IllegalStateException("anchor end reference pos must be either anchor start or anchor end")
            }
        }

        // use the mate mapped position to see if we are relevant
        // this function is a little bit more complicated due to that
        // Essentially, we want to find the following situation, where
        // we have an unmapped read and we want the mate to be up stream for V, and downstream for J
        //
        // >------------V--------D----------J---------> positive strand
        // ======>  <=====              ====>  <=====
        //  mate     this               this    mate
        //
        // or
        //
        // <------------J--------D----------V---------< negative strand
        // ======>  <=====              ====>  <=====
        //  mate     this                this    mate
        //
        fun isUnamppedReadRelevantToAnchorLoc(
            read: SAMRecord, anchorLocation: VJAnchorGenomeLocation, maxFragmentLength: Int): Boolean
        {
            if (!read.readPairedFlag || read.mateUnmappedFlag)
                return false

            if (anchorLocation.chromosome != read.mateReferenceName)
                return false

            val mateMappedStrand: Strand = if (read.mateNegativeStrandFlag) Strand.REVERSE else Strand.FORWARD
            val mateMappedStart: Int = read.mateAlignmentStart

            // we cannot easily get the mate alignment end, so we just infer it for now
            val inferredMateMappedEnd: Int = mateMappedStart + read.readLength

            if (anchorLocation.vj == VJ.V && anchorLocation.strand == Strand.FORWARD ||
                anchorLocation.vj == VJ.J && anchorLocation.strand == Strand.REVERSE)
            {
                // here we want the anchor to have higher coord than mate read, so the mate read must be on forward strand
                if (mateMappedStrand != Strand.FORWARD)
                    return false

                if (mateMappedStart > anchorLocation.end)
                    return false

                // next check that we are not too far away
                return (anchorLocation.start - inferredMateMappedEnd) < maxFragmentLength
            }
            else
            {
                // here we want the anchor to have lower coord than mate read, so the read must be on reverse strand
                if (mateMappedStrand != Strand.REVERSE)
                    return false

                if (inferredMateMappedEnd < anchorLocation.start)
                    return false

                // next check that we are not too far away
                return (mateMappedStart - anchorLocation.end) < maxFragmentLength
            }
        }

        // use the mate mapped position to see if we are relevant to a constant region
        // we have an unmapped read and we want the mate to be mapped to C and pointing up stream
        //
        // >-----V--------D----------J----------C------> positive strand
        //                        ======>      <=====
        //                         this         mate
        //
        // or
        //
        // <-----C-------J--------D----------V---------< negative strand
        //   ======>    <=====
        //    mate       this
        //
        fun isUnamppedReadRelevantToConstantRegion(
            read: SAMRecord, constantRegionLocation: GenomicLocation, maxFragmentLength: Int): Boolean
        {
            if (!read.readPairedFlag || read.mateUnmappedFlag)
                return false

            if (constantRegionLocation.chromosome != read.mateReferenceName)
                return false

            val mateMappedStart: Int = read.mateAlignmentStart

            // we cannot easily get the mate alignment end, so we just infer it for now
            val inferredMateMappedEnd: Int = mateMappedStart + read.readLength

            // too far away
            if ((mateMappedStart - constantRegionLocation.posEnd) >= maxFragmentLength ||
                (constantRegionLocation.posStart - inferredMateMappedEnd) >= maxFragmentLength)
            {
                return false
            }

            val mateMappedStrand: Strand = if (read.mateNegativeStrandFlag) Strand.REVERSE else Strand.FORWARD

            // see picture above, the strand must not be equal
            return constantRegionLocation.strand != mateMappedStrand
        }
    }
}