package com.hartwig.hmftools.esvee.alignment;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;

public class AlignData
{
    public final ChrBaseRegion RefLocation;
    public final int SequenceStart;
    public final int SequenceEnd;
    public final int MapQual;
    public final int NMatches;
    public final int Score;
    public final String Cigar;
    public final String LocationInfo;
    public final String MdTag;

    public AlignData(
            final ChrBaseRegion refLocation, final int sequenceStart, final int sequenceEnd, final int mapQual,
            final int nMatches, final int score, final String cigar, final String locationInfo, final String mdTag)
    {
        RefLocation = refLocation;
        SequenceStart = sequenceStart;
        SequenceEnd = sequenceEnd;
        MapQual = mapQual;
        NMatches = nMatches;
        Score = score;
        Cigar = cigar;
        LocationInfo = locationInfo;
        MdTag = mdTag;
    }

    public static AlignData from(final BwaMemAlignment alignment, final RefGenomeVersion refGenomeVersion)
    {
        int chrIndex = alignment.getRefId();

        if(chrIndex < 0 || chrIndex >= HumanChromosome.values().length)
            return null;

        String chromosome = refGenomeVersion.versionedChromosome(HumanChromosome.values()[chrIndex].toString());

        return new AlignData(
                new ChrBaseRegion(chromosome, alignment.getRefStart(), alignment.getRefEnd()),
                alignment.getSeqStart(), alignment.getSeqEnd(), alignment.getMapQual(), alignment.getAlignerScore(),
                alignment.getNMismatches(), alignment.getCigar(), alignment.getXATag(), alignment.getMDTag());
    }

    public String toString()
    {
        return String.format("%s %s seq(%d-%d) score(%d) mq(%d)",
                RefLocation, Cigar, SequenceStart, SequenceEnd, Score, MapQual);
    }
}
