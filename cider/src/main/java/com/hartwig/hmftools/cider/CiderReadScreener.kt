package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.CiderUtils.getAdapterDnaTrim
import com.hartwig.hmftools.cider.CiderUtils.setAdapterDnaTrim
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
import kotlin.math.abs
import kotlin.math.min

class CiderReadScreener(// collect the reads and sort by types
    private val mCiderGeneDatastore: ICiderGeneDatastore,
    private val mAnchorBlosumSearcher: IAnchorBlosumSearcher,
    private val mMinReadCdr3Overlap: Int,
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

    private val mProcessedReadRecords = Collections.newSetFromMap(ConcurrentHashMap<ReadRecordKey, Boolean>())
    private val mVjReadCandidates = Collections.synchronizedList(ArrayList<VJReadCandidate>())
    private val mAllMatchedReads = Collections.synchronizedList(ArrayList<SAMRecord>())
    private val mMatchedReadNames = Collections.newSetFromMap(ConcurrentHashMap<String, Boolean>())
    private val mUnmatchedMappedReads = Collections.synchronizedList(ArrayList<SAMRecord>())
    private val mVjMateCandidates = ArrayList<VJReadCandidate>()

    val vjReadCandidates: List<VJReadCandidate>
        get() = mVjReadCandidates
    val allMatchedReads: List<SAMRecord>
        get() = mAllMatchedReads
    val vjMateCandidates: List<VJReadCandidate>
        get() = mVjMateCandidates

    // Note: this function is called from multiple threads
    fun asyncProcessSamRecord(samRecord: SAMRecord)
    {
        val readRecordKey = ReadRecordKey(samRecord.readName, SamRecordUtils.firstInPair(samRecord),
                                            samRecord.referenceIndex, samRecord.alignmentStart)

        if (!mProcessedReadRecords.add(readRecordKey))
        {
            //sLogger.trace("read record {} already processed", readRecordKey);
            return
        }

        // Calculate the trimming required to remove any sequencing adapter bases left on the read (for short fragments).
        // The trimming amount is stored into the read so we can retrieve it later when needed.
        setAdapterDnaTrim(samRecord)

        tryMatchRead(samRecord)
    }

    // Collect mates of candidate reads, which can be used to extend the consensus sequence.
    fun processCandidateMates()
    {
        // For each unmatched read, try to align it to the anchor(s) which matched to its fragment from mates.
        // Only the alignment matching method is used because other methods would've matched on the first pass.

        val anchorLocationsByFragment = mVjReadCandidates
            .groupBy { it.read.readName }
            .mapValues {
                (_, candidates) -> candidates.flatMap {
                    candidate -> candidate.vjAnchorTemplates.mapNotNull {
                        template -> template.anchorLocation?.let {
                            location -> VJAnchorGenomeLocation(template.type, location) } } } }
        for (read in mUnmatchedMappedReads)
        {
            // Since we are matching with tolerance up to the fragment length, we need to choose the closest match.
            // Some genes are within a few hundred bases of each other, so the read could match multiple.
            val fragmentAnchorLocations = anchorLocationsByFragment[read.readName] ?: continue
            val mapped = GenomeRegions.create(read.referenceName, read.alignmentStart, read.alignmentEnd)
            val readCandidates = fragmentAnchorLocations.mapNotNull { anchorLocation ->
                tryMatchToAnchorLocation(read, mapped, anchorLocation, mMaxFragmentLength)
                    ?.let { readCandidate ->
                        val distanceFromAnchor = min(abs(readCandidate.anchorOffsetStart), abs(readCandidate.anchorOffsetEnd))
                        Pair(readCandidate, distanceFromAnchor)
                    }
            }
            readCandidates.minByOrNull { it.second }?.let { addVjReadCandidate(it.first, true) }
        }
    }

    private fun tryMatchRead(samRecord: SAMRecord)
    {
        val fragmentMapped = samRecord.referenceIndex != -1
        val readMapped = !samRecord.readUnmappedFlag

        val readCandidate = if (fragmentMapped)
            (if (readMapped) tryMatchByAlignment(samRecord)
            // if read is not mapped then the mate is mapped to this region instead, we use the mate mapped
            // region to decide which locus (IGH, TRA etc) we need to search for
            // if we do not do this then we will very likely match to wrong type
            else tryMatchUnmappedReadByBlosum(samRecord))
        // if both read and its mate are unmapped then we just try match by blosum
        else tryMatchByBlosum(samRecord)

        if (readCandidate != null)
        {
            mAllMatchedReads.add(samRecord)
            mMatchedReadNames.add(samRecord.readName)
            addVjReadCandidate(readCandidate, false)
        }
        else if (readMapped)
        {
            // If the mate of this read is relevant, then later we may include this read too.
            // Only keep mapped reads because we will need the aligned position.
            mUnmatchedMappedReads.add(samRecord)
        }
    }

    private fun tryMatchByAlignment(samRecord: SAMRecord): VJReadCandidate?
    {
        val mapped = GenomeRegions.create(
            samRecord.referenceName,
            samRecord.alignmentStart, samRecord.alignmentEnd)

        // first step we see if the read overlaps with any anchor location
        tryMatchToAnchorLocationExactly(samRecord, mapped)
            ?.let { return it }

        val leftSoftClip = CigarUtils.leftSoftClipLength(samRecord)
        val rightSoftClip = CigarUtils.rightSoftClipLength(samRecord)
        if (leftSoftClip != 0 || rightSoftClip != 0)
        {
            for (igTcrConstantRegion in mCiderGeneDatastore.getIgConstantDiversityRegions())
            {
                // now try to match around location of constant regions
                tryMatchFromConstantRegion(samRecord, mapped, igTcrConstantRegion)
                    ?.let { return it }
            }
        }
        return null
    }

    private fun tryMatchToAnchorLocationExactly(samRecord: SAMRecord, mapped: GenomeRegion): VJReadCandidate?
    {
        for (anchorLocation in mCiderGeneDatastore.getVjAnchorGeneLocations())
        {
            tryMatchToAnchorLocation(samRecord, mapped, anchorLocation, 0)
                ?.let { return it }
        }
        return null
    }

    internal fun tryMatchToAnchorLocation(
        samRecord: SAMRecord, mapped: GenomeRegion,
        anchorLocation: VJAnchorGenomeLocation,
        tolerance: Int
    ): VJReadCandidate?
    {
        if (!isMappedToAnchorLocation(mapped, anchorLocation, tolerance))
        {
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

        val readTrim = getAdapterDnaTrim(samRecord)

        // read:   |---------------|    |---------------------|
        // anchor:    |----|                   |----|
        // CDR3:           |-------------------|
        // check:          <------->    <------>
        val isNearbyAnchor =
                (anchorLocation.anchorBoundarySide() == 1 && samRecord.readLength - readTrim.second - readAnchorEnd >= mMinReadCdr3Overlap - tolerance) ||
                (anchorLocation.anchorBoundarySide() == -1 && readAnchorStart - readTrim.first >= mMinReadCdr3Overlap - tolerance)

        if (isNearbyAnchor)
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
            val genes = mCiderGeneDatastore.getByGeneLocation(anchorLocation.genomeLocation)

            if (genes.isNotEmpty())
            {
                return createVjReadCandidate(samRecord, genes, MatchMethod.ALIGN,
                    anchorLocation.strand === Strand.REVERSE,
                    readAnchorStart, readAnchorEnd, anchorLocation.genomeLocation)
            }
        }
        return null
    }

    // This read is unmapped and mate is mapped. We use the mate mapping to find out
    // which Ig/TCR locus this read is near. And only blosum search for anchor with
    // this locus. This is needed to make sure we find the correct locus
    private fun tryMatchUnmappedReadByBlosum(read: SAMRecord): VJReadCandidate?
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
                return null
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
            return null
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
            if (anchorBlosumMatch != null && anchorBlosumMatch.similarityScore > 0 && anchorBlosumMatch.templateGenes.isNotEmpty())
            {
                val readCandidate = createVjReadCandidate(
                    read,
                    anchorBlosumMatch.templateGenes,
                    MatchMethod.BLOSUM,
                    reverseRead,
                    anchorBlosumMatch.anchorStart,
                    anchorBlosumMatch.anchorEnd,
                    null
                )
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
                return readCandidate
            }
        }

        return null
    }

    private fun tryMatchByBlosum(samRecord: SAMRecord): VJReadCandidate?
    {
        // in case someone changes this
        assert(Strand.entries.size == 2)

        for (strand in Strand.entries)
        {
            var readString = samRecord.readString
            if (strand == Strand.REVERSE) readString = SequenceUtil.reverseComplement(readString)

            val anchorBlosumMatch: AnchorBlosumMatch? =
                mAnchorBlosumSearcher.searchForAnchor(readString, IAnchorBlosumSearcher.Mode.DISALLOW_NEG_SIMILARITY)

            if (anchorBlosumMatch != null && anchorBlosumMatch.similarityScore > 0 && anchorBlosumMatch.templateGenes.isNotEmpty())
            {
                return createVjReadCandidate(
                    samRecord,
                    anchorBlosumMatch.templateGenes,
                    MatchMethod.BLOSUM,
                    strand == Strand.REVERSE,
                    anchorBlosumMatch.anchorStart,
                    anchorBlosumMatch.anchorEnd,
                    null
                )
            }
        }
        return null
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
            if (anchorBlosumMatch != null && anchorBlosumMatch.similarityScore > 0 && anchorBlosumMatch.templateGenes.isNotEmpty())
            {
                sLogger.trace("read({}) matched from constant region({}) using blosum", samRecord, igTcrConstantDiversityRegion)

                return createVjReadCandidate(
                    samRecord,
                    anchorBlosumMatch.templateGenes,
                    MatchMethod.BLOSUM,
                    strand == Strand.REVERSE,
                    anchorBlosumMatch.anchorStart,
                    anchorBlosumMatch.anchorEnd,
                    null
                )
            }
        }
        return null
    }

    private fun createVjReadCandidate(
        samRecord: SAMRecord,
        vjAnchorTemplates: List<VJAnchorTemplate>,
        templateMatchMethod: MatchMethod,
        useRevComp: Boolean,
        readAnchorStart: Int,
        readAnchorEnd: Int,
        templateLocation: GenomicLocation?)
    : VJReadCandidate
    {
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
            throw RuntimeException("blosum match with -ve similarity score: $readMatch")
        }

        return readMatch
    }

    private fun addVjReadCandidate(candidate: VJReadCandidate, isMate: Boolean)
    {
        if (isMate)
        {
            mVjMateCandidates.add(candidate)
        }
        else
        {
            mVjReadCandidates.add(candidate)
        }
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(CiderReadScreener::class.java)

        fun isMappedToAnchorLocation(mapped: GenomeRegion, anchorLocation: VJAnchorGenomeLocation, tolerance: Int): Boolean
        {
            if (anchorLocation.chromosome != mapped.chromosome()) return false
            return anchorLocation.start <= mapped.end() + tolerance && mapped.start() - tolerance <= anchorLocation.end
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