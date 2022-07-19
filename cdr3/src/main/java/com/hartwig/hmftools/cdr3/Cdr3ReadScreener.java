package com.hartwig.hmftools.cdr3;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import javax.annotation.Nullable;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.samtools.SamRecordUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.SequenceUtil;

public class Cdr3ReadScreener
{
    private static final Logger sLogger = LogManager.getLogger(Cdr3ReadScreener.class);

    // collect the reads and sort by types
    private VJGeneStore mVJGeneStore;
    int mMaxAnchorAlignDistance;

    private Multimap<ReadKey, VJReadCandidate> mCdr3ReadMatchMap = ArrayListMultimap.create();

    // want a map to make sure we do not process same record twice
    private final Set<VJReadRecordKey> mProcessedReadRecords = new HashSet<>();

    private final List<SAMRecord> mAllMatchedReads = new ArrayList<>();

    public Cdr3ReadScreener(VJGeneStore vjGeneStore, int maxAnchorAlignDistance)
    {
        mVJGeneStore = vjGeneStore;
        mMaxAnchorAlignDistance = maxAnchorAlignDistance;
    }

    public Multimap<ReadKey, VJReadCandidate> getVJReadCandidates()
    {
        return mCdr3ReadMatchMap;
    }

    public List<SAMRecord> getAllMatchedReads()
    {
        return mAllMatchedReads;
    }

    public void processSamRecord(SAMRecord samRecord)
    {
        // see if we already processed this read
        // note that we cannot use unmapped flag, we must check the reference index against -1
        // this is cause in many cases read record will use
        GenomeRegion mapped = samRecord.getReferenceIndex() == -1 ? null : GenomeRegions.create(samRecord.getReferenceName(),
                samRecord.getAlignmentStart(), samRecord.getAlignmentEnd());

        VJReadRecordKey readRecordKey = new VJReadRecordKey(samRecord.getReadName(), samRecord.getFirstOfPairFlag(), mapped);

        if (!mProcessedReadRecords.add(readRecordKey))
        {
            sLogger.trace("read record {} already processed", readRecordKey);
            return;
        }

        Set<VJGeneType> matchedGeneTypes = new HashSet<>();

        boolean matchFound = false;

        // first we try to match by genome region
        if (mapped != null)
        {
            matchFound |= tryMatchByAlignment(samRecord, mapped, matchedGeneTypes);
        }

        // next match by direct compare, we always have to do both, reason is that if this is a
        // CDR3 gene then a read could match V, D, or J
        matchFound |= tryMatchByExact(samRecord, matchedGeneTypes);

        if (matchFound)
        {
            mAllMatchedReads.add(samRecord);
        }
        else
        {
            sLogger.trace("read record {} does not match any VJ", samRecord);
        }
    }

    private boolean tryMatchByAlignment(final SAMRecord samRecord, final GenomeRegion mapped, Set<VJGeneType> matchedGenes)
    {
        boolean anchorFound = false;

        for (GeneLocation anchorLocation : mVJGeneStore.getVJGenomeRegions())
        {
            // see if the anchor location is mapped around here
            if ((anchorLocation.start() - mMaxAnchorAlignDistance <= mapped.end() &&
                anchorLocation.end() + mMaxAnchorAlignDistance >= mapped.start()) &&
                anchorLocation.chromosome().equals(mapped.chromosome()))
            {
                int anchorLength = anchorLocation.baseLength();

                if (anchorLength != 30)
                {
                    throw new RuntimeException("unexpected anchor length");
                }

                // print out the read id and sequence, but first we need to know which gene it came from
                // print the anchor sequence
                // use 0 base positions, easier for strings
                int readAnchorStart = extrapolateReadOffsetAtRefPosition(samRecord, anchorLocation.start(), mapped);
                int readAnchorEnd = extrapolateReadOffsetAtRefPosition(samRecord, anchorLocation.end(), mapped);

                if (readAnchorStart == -1 || readAnchorEnd == -1)
                {
                    // anchor pos not here
                    continue;
                }
                readAnchorEnd += 1;

                if (readAnchorStart >= 0 && readAnchorEnd >= 0 && readAnchorEnd <= samRecord.getReadLength())
                {
                    if (anchorLocation.getStrand() == Strand.REVERSE)
                    {
                        // say if read length is 2, and read anchor start is 1, read anchor end is 2
                        // reverse it we want start = 1, end = 2
                        int revStart = samRecord.getReadLength() - readAnchorEnd;
                        readAnchorEnd = samRecord.getReadLength() - readAnchorStart;
                        readAnchorStart = revStart;
                    }

                    // we want to make sure same gene is not included twice
                    List<VJGene> genes = mVJGeneStore.getByAnchorGeneLocation(anchorLocation).stream()
                            .filter(o -> !matchedGenes.contains(o))
                            .collect(Collectors.toList());

                    if (!genes.isEmpty())
                    {
                        matchedGenes.addAll(genes.stream().map(VJGene::getType).collect(Collectors.toList()));

                        addCdr3ReadMatch(samRecord,
                                genes,
                                VJReadCandidate.AnchorMatchType.ALIGN,
                                anchorLocation.getStrand() == Strand.REVERSE,
                                readAnchorStart,
                                readAnchorEnd,
                                anchorLocation);
                        anchorFound = true;
                    }
                }
            }
        }
        return anchorFound;
    }

    public boolean tryMatchByExact(SAMRecord samRecord, Set<VJGeneType> matchedGenes)
    {
        String seq = samRecord.getReadString();
        String seqRevComp = SequenceUtil.reverseComplement(seq);
        boolean useRevComp = false;
        boolean matchFound = false;

        // for each record we need to find the anchor sequence and see if it matches
        for (String anchorSeq : mVJGeneStore.getAnchorSequenceSet())
        {
            int anchorIndex = seq.indexOf(anchorSeq);

            if (anchorIndex != -1)
            {
                useRevComp = false;
            }
            else if ((anchorIndex = seqRevComp.indexOf(anchorSeq)) != -1)
            {
                useRevComp = true;
            }

            if (anchorIndex != -1)
            {
                // the anchor sequence aligns with frame so we need to make that adjustment
                String geneSeq = seq;

                if (useRevComp)
                {
                    geneSeq = SequenceUtil.reverseComplement(geneSeq);
                }

                // want to make sure same gene is not included twice
                List<VJGene> genes = mVJGeneStore.getByAnchorSequence(anchorSeq).stream()
                        .filter(o -> !matchedGenes.contains(o.getType()))
                        .collect(Collectors.toList());

                if (!genes.isEmpty())
                {
                    addCdr3ReadMatch(samRecord,
                            genes,
                            VJReadCandidate.AnchorMatchType.EXACT,
                            useRevComp,
                            anchorIndex,
                            anchorIndex + anchorSeq.length(),
                            null);

                    matchFound = true;

                /*sLogger.info("anchor found but does not overlap any gene: read: id={} isSupp={}, anchor seq={}," +
                                " mapped={}:{}-{}, gene seq={}, AA seq={}",
                        samRecord.getReadName(), samRecord.isSecondaryOrSupplementary(), anchorSeq,
                        samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(),
                        geneSeq, aaSeq);*/
                }
            }
        }

        return matchFound;
    }

    public void addCdr3ReadMatch(SAMRecord samRecord, List<VJGene> VJGenes,
            VJReadCandidate.AnchorMatchType anchorMatchType, boolean useRevComp,
            int readAnchorStart, int readAnchorEnd, @Nullable GeneLocation anchorLocation)
    {
        if (VJGenes.isEmpty())
            return;

        String seq = samRecord.getReadString();

        if (useRevComp)
            seq = SequenceUtil.reverseComplement(seq);

        String anchorSeq = seq.substring(readAnchorStart, readAnchorEnd);

        // find out the imgt gene type. They should be the same type
        VJGeneType geneType = VJGenes.stream().findFirst().get().getType();

        // check to make sure all the same
        if (VJGenes.stream().anyMatch(o -> o.getType() != geneType))
        {
            sLogger.error("multiple gene types found in same match: {}", VJGenes);
            throw new RuntimeException("multiple gene types found in same match");
        }

        // since we don't actually know whether the aligned part is the anchor sequence, we have to use
        // the soft clip that we think make sense
        int leftSoftClip = SamRecordUtils.leftSoftClip(samRecord);
        int rightSoftClip = SamRecordUtils.rightSoftClip(samRecord);

        VJReadCandidate readMatch = new VJReadCandidate(samRecord, VJGenes, geneType, anchorMatchType, useRevComp,
                readAnchorStart, readAnchorEnd, anchorLocation, leftSoftClip, rightSoftClip);

        mCdr3ReadMatchMap.put(new ReadKey(samRecord.getReadName(), samRecord.getFirstOfPairFlag()), readMatch);

        List<String> geneNames = VJGenes.stream().map(o -> o.getName()).distinct().collect(Collectors.toList());

        sLogger.info("genes: {}, read: {}, match type={} anchor seq={}, soft clipped={};{}, gene seq={}, AA seq={}",
                geneNames, samRecord, anchorMatchType, anchorSeq,
                leftSoftClip, rightSoftClip,
                readMatch.getPotentialCdr3Dna(), readMatch.getPotentialCdr3AA());
    }

    static int extrapolateReadOffsetAtRefPosition(SAMRecord samRecord, int refPos, final GenomeRegion mapped)
    {
        // print out the read id and sequence, but first we need to know which gene it came from
        // print the anchor sequence
        // use 0 base positions, easier for strings
        int readOffset = samRecord.getReadPositionAtReferencePosition(refPos, false) - 1;
        if (readOffset >= 0)
        {
            return readOffset;
        }

        // cannot find it, we try to extrapolate
        if (mapped.end() < refPos)
        {
            // the mapped location is below the anchor location
            // we want to see if this read bases extends all the way to the anchor start
            int mappedEndOffset = samRecord.getReadPositionAtReferencePosition(mapped.end(), false) - 1;
            readOffset = mappedEndOffset + refPos - mapped.end();
        }
        else if (mapped.start() > refPos)
        {
            // the mapped location is above the anchor location
            // we want to see if this read bases extends all the way to the anchor end
            int mappedStartOffset = samRecord.getReadPositionAtReferencePosition(mapped.start(), false) - 1;
            readOffset = mappedStartOffset + refPos - mapped.start();
        }
        else
        {
            // shouldn't get here
            assert(false);
        }
        if (readOffset < 0 || readOffset >= samRecord.getReadLength())
        {
            return -1;
        }

        return readOffset;
    }
}
