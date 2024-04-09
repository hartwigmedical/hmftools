package com.hartwig.hmftools.esvee.alignment;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import static htsjdk.samtools.SAMFlag.READ_REVERSE_STRAND;

import com.hartwig.hmftools.common.bam.SamRecordUtils;
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
    public final int Flags;
    public final String Cigar;
    public final String XaTag;
    public final String MdTag;

    public AlignData(
            final ChrBaseRegion refLocation, final int sequenceStart, final int sequenceEnd, final int mapQual,
            final int score, final int flags, final String cigar, final int nMatches, final String xaTag, final String mdTag)
    {
        RefLocation = refLocation;
        SequenceStart = sequenceStart;
        SequenceEnd = sequenceEnd;
        MapQual = mapQual;
        NMatches = nMatches;
        Score = score;
        Flags = flags;

        Cigar = cigar;
        XaTag = xaTag;
        MdTag = mdTag;
    }

    public byte orientation()
    {
        return SamRecordUtils.isFlagSet(Flags, READ_REVERSE_STRAND) ? NEG_ORIENT : POS_ORIENT;
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
                alignment.getSamFlag(), alignment.getCigar(), alignment.getNMismatches(), alignment.getXATag(), alignment.getMDTag());
    }

    public String toString()
    {
        return String.format("%s %s seq(%d-%d) score(%d) mq(%d) flags(%d)",
                RefLocation, Cigar, SequenceStart, SequenceEnd, Score, MapQual, Flags);
    }
}
