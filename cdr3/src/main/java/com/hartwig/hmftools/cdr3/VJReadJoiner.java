package com.hartwig.hmftools.cdr3;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.codon.Codons;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// try to see if we can join up V and J reads together
public class VJReadJoiner
{
    public static final Logger sLogger = LogManager.getLogger(VJReadJoiner.class);

    private Multimap<ReadKey, VJReadCandidate> mCdr3ReadMatchMap;

    public VJReadJoiner(Multimap<ReadKey, VJReadCandidate> cdr3ReadCandidates)
    {
        // we need to sort all the reads into V and J

        // for now have to go through all the matches and tease out the reads

        Map<ReadKey, VJReadCandidate> ighvReadCandidates = findReadMatchOfType(cdr3ReadCandidates, VJGeneType.IGHV);
        Map<ReadKey, VJReadCandidate> ighjReadCandidates = findReadMatchOfType(cdr3ReadCandidates, VJGeneType.IGHJ);

        // find reads that contain both V and J, they should form a VDJ
        Set<ReadKey> readKeys = new HashSet<>(ighvReadCandidates.keySet());
        readKeys.retainAll(ighjReadCandidates.keySet());

        sLogger.info("{} reads with both V and J match found", readKeys.size());

        for (ReadKey readKey : readKeys)
        {
            // each of these reads has both a V and a J section, they are VDJ candidates, for now just log them
            // Must satisfy 2 rules:
            // 1. the V and J have to be in same direction
            // 2. V has to come before J
            // 3. the D section should be of sufficient length

            VJReadCandidate vMatch = ighvReadCandidates.get(readKey);
            VJReadCandidate jMatch = ighjReadCandidates.get(readKey);

            combineVJMatch(readKey, vMatch, jMatch);
        }
    }

    public static Cdr3ReadVJMatch combineVJMatch(ReadKey readKey, @Nullable VJReadCandidate vMatch, @Nullable VJReadCandidate jMatch)
    {
        if (vMatch == null && jMatch == null)
            return null;

        if (vMatch != null && jMatch == null)
        {
            return new Cdr3ReadVJMatch(vMatch, null, vMatch.getReadSequence(), 0, -1);
        }

        if (vMatch == null && jMatch != null)
        {
            return new Cdr3ReadVJMatch(null, jMatch, jMatch.getReadSequence(), -1, 0);
        }

        // each of these reads has both a V and a J section, they are VDJ candidates, for now just log them
        // Must satisfy 2 rules:
        // 1. the V and J have to be in same direction
        // 2. V has to come before J
        // 3. the D section should be of sufficient length

        sLogger.info("read {} vSeq: {} jSeq: {} vRead: {} {} jRead: {} {}",
                readKey, vMatch.getReadSequence(), jMatch.getReadSequence(),
                vMatch.getRead(), vMatch.getRead().getCigarString(),
                jMatch.getRead(), jMatch.getRead().getCigarString());

        // they should be of the same read already
        // first we have to establish the direction of the read
        if ((vMatch.getRead().getReadNegativeStrandFlag() == jMatch.getRead().getReadNegativeStrandFlag()) !=
                (vMatch.getUseReverseComplement() == jMatch.getUseReverseComplement()))
        {
            sLogger.debug("skipping read {} v and j direction mismatch", readKey);
            return null;
        }

        // if one is supplementary the other is not, we have to adjust the offset and find the "true" sequence
        String vReadSeq = vMatch.getReadSequence();
        String jReadSeq = jMatch.getReadSequence();

        String seq = vReadSeq;
        int vSeqOffset = 0;
        int jSeqOffset = 0;

        if (vReadSeq.length() < jReadSeq.length())
        {
            seq = jReadSeq;
            vSeqOffset = seq.indexOf(vReadSeq);
        }
        else if (vReadSeq.length() > jReadSeq.length())
        {
            jSeqOffset = seq.indexOf(jReadSeq);
        }

        if (vSeqOffset == -1 || jSeqOffset == -1)
        {
            sLogger.warn("Cannot process read {}, unable to determine v or j offset. vSeq: {} jSeq: {} vRead: {} {} jRead: {} {}",
                    readKey, vMatch.getReadSequence(), jMatch.getReadSequence(),
                    vMatch.getRead(), vMatch.getRead().getCigarString(),
                    jMatch.getRead(), jMatch.getRead().getCigarString());
            return null;
        }

        if (vMatch.getAnchorOffsetEnd() + vSeqOffset >= jMatch.getAnchorOffsetStart() + jSeqOffset)
        {
            sLogger.debug("skipping read {} v and j anchor location incorrect", readKey);
            return null;
        }


        // now read out the sequence

        String vdjSeq = seq.substring(vMatch.getAnchorOffsetStart() + vSeqOffset, jMatch.getAnchorOffsetEnd() + jSeqOffset);

        sLogger.info("read {} vdj DNA seq: {} vdj amino acid seq: {}", readKey, vdjSeq, Codons.aminoAcidFromBases(vdjSeq));
        //sLogger.info("read {} vSeq: {} jSeq: {}", readKey, vMatch.getReadSequence(), jMatch.getReadSequence());

        return new Cdr3ReadVJMatch(vMatch, jMatch, seq, vSeqOffset, jSeqOffset);
    }

    private static Map<ReadKey, VJReadCandidate> findReadMatchOfType(
            Multimap<ReadKey, VJReadCandidate> cdr3ReadMatchMap,
            VJGeneType geneType)
    {
        Map<ReadKey, VJReadCandidate> readCandidates = new HashMap<>();

        for (var e : cdr3ReadMatchMap.entries())
        {
            ReadKey key = e.getKey();
            VJReadCandidate match = e.getValue();
            if (match.getVjGeneType() == geneType)
            {
                VJReadCandidate existingMatch = readCandidates.get(key);
                if (existingMatch != null)
                {
                    //sLogger.warn("Conflicting match of type {} found: {} vs {}", geneType, existingMatch, match);
                }

                readCandidates.put(e.getKey(), match);
            }
        }

        return readCandidates;
    }
}
