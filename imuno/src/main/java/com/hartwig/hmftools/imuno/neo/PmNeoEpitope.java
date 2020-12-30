package com.hartwig.hmftools.imuno.neo;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.switchStream;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.ALL_TRANS_BASES;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.getDownstreamCodingBases;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.getUpstreamCodingBases;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.setTranscriptCodingData;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.setTranscriptContext;

import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class PmNeoEpitope extends NeoEpitope
{
    private final PointMutationData mPointMutation;
    private final int mIndelBaseDiff;

    public PmNeoEpitope(final PointMutationData pointMutation)
    {
        super();
        mPointMutation = pointMutation;
        mIndelBaseDiff = mPointMutation.Alt.length() - mPointMutation.Ref.length();
    }

    public byte orientation(int fs)
    {
        return (fs == FS_UP) == (TransData[fs].Strand == POS_STRAND) ? POS_ORIENT : NEG_ORIENT;
    }

    public int position(int stream)
    {
        // position for a point mutation is defined on the upstream as the base of the mutation
        // for the downstream it is one base before for +ve strand and 1 after for -ve strand
        int pmPosition = mPointMutation.Position;

        if(mIndelBaseDiff >= 0)
        {
            if(stream == FS_UP)
                return pmPosition;
            else
                return TransData[FS_UP].Strand == POS_STRAND ? pmPosition + 1 : pmPosition - 1;
        }

        // handle the DEL scenario - shift the downstream position by the delete length
        if(stream == FS_UP)
        {
            return pmPosition;
        }
        else
        {
            if(TransData[FS_UP].Strand == POS_STRAND)
                return pmPosition + abs(mIndelBaseDiff) + 1;
            else
                return pmPosition - 1;
        }
    }

    public String chromosome(int stream)
    {
        return mPointMutation.Chromosome;
    }

    public String geneName(int stream) { return mPointMutation.Gene; }

    public String variantType()
    {
        // missense, inframe insertion, inframe deletion, frameshift,
        if(mPointMutation.Effect == MISSENSE || mIndelBaseDiff == 0)
            return MISSENSE.toString();

        if((mIndelBaseDiff % 3) == 0)
        {
            if(mIndelBaseDiff > 0)
                return "INFRAME_INSERTION";
            else
                return "INFRAME_DELETION";
        }

        return "FRAMESHIFT";
    }

    public String variantInfo()
    {
        return String.format("%s:%d:%s:%s", mPointMutation.Chromosome,  mPointMutation.Position, mPointMutation.Ref, mPointMutation.Alt);
    }

    public double copyNumber() { return mPointMutation.CopyNumber; }

    public void setTranscriptData(final TranscriptData upTransData, final TranscriptData downTransData)
    {
        final TranscriptData transData = upTransData;

        TransData[FS_UP] = TransData[FS_DOWN] = transData;

        int indelBaseDiff = mPointMutation.Alt.length() - mPointMutation.Ref.length();
        int insertedBaseLength = max(indelBaseDiff, 0);

        // set the data for the lower part of the mutation
        int lowerStream = transData.Strand == POS_STRAND ? FS_UP : FS_DOWN;

        setTranscriptContext(this, transData, position(lowerStream), lowerStream);

        int insSeqLength = lowerStream == FS_UP ? insertedBaseLength : 0;

        setTranscriptCodingData(this, transData, position(lowerStream), insSeqLength, lowerStream);

        int upperStream = switchStream(lowerStream);

        // for DELs, set the downstream data as well since it can cross exon-boundaries and/or affect coding bases
        setTranscriptContext(this, transData, position(upperStream), upperStream);

        insSeqLength = upperStream == FS_UP ? insertedBaseLength : 0;
        setTranscriptCodingData(this, transData, position(upperStream), insSeqLength, upperStream);
    }

    public int getDownstreamPhaseOffset()
    {
        int upOpenBases = getUpstreamOpenCodonBases();

        if(upOpenBases == 0)
        {
            if(mPointMutation.Effect == MISSENSE)
                return 3;
            else
                return 0;
        }

        return 3 - upOpenBases;
    }

    public void extractCodingBases(final RefGenomeInterface refGenome, int requiredAminoAcids)
    {
        // get the number of bases from the mutation or junction as required by the required amino acid count
        // the downstream bases start directly after the mutation (ie after any insert, or deleted bases)

        // if the upstream phase is not 2 (ie the end of a codon) then additionally extract enough bases to cover the required amino
        // acids plus the bases of the partial codon - the downstream bases will provide the remainder to complete the codon

        // phasing already takes any inserted bases (ie from an INDEL) into account, so just adjust for the number of bases required
        int upPhaseOffset = getUpstreamOpenCodonBases();
        int downPhaseOffset = getDownstreamPhaseOffset();

        int upRequiredBases = requiredAminoAcids * 3 + upPhaseOffset;

        String upCodingBases = "";
        String downCodingBases = "";

        upRequiredBases -= mPointMutation.Alt.length();

        int deletedBases = mIndelBaseDiff < 0 ? abs(mIndelBaseDiff) : 0;

        int upPosition;
        int downPosition;

        if(TransData[FS_UP].Strand == POS_STRAND)
        {
            upPosition = mPointMutation.Position - 1;
            downPosition = mPointMutation.Position + deletedBases + 1;
        }
        else
        {
            upPosition = mPointMutation.Position + deletedBases + 1;
            downPosition = mPointMutation.Position - 1;
        }

        upCodingBases = getUpstreamCodingBases(
                refGenome, TransData[FS_UP], chromosome(FS_UP), upPosition, orientation(FS_UP), upRequiredBases);

        // the ref bases was excluded and is now replaced by the alt
        if(TransData[FS_UP].Strand == POS_STRAND)
        {
            upCodingBases += mPointMutation.Alt;
        }
        else
        {
            upCodingBases = mPointMutation.Alt + upCodingBases;
        }

        boolean canStartInExon = true; // assumed true for now RegionType[FS_UPSTREAM] == TranscriptRegionType.EXONIC;
        int downRequiredBases = phaseMatched() ? requiredAminoAcids * 3 + downPhaseOffset : ALL_TRANS_BASES;

        downCodingBases = getDownstreamCodingBases(
                refGenome, TransData[FS_DOWN], chromosome(FS_DOWN), downPosition, orientation(FS_DOWN),
                downRequiredBases, canStartInExon, false);

        CodingBases[FS_UP] = upCodingBases;
        CodingBases[FS_DOWN] = downCodingBases;
    }

    public String toString()
    {
        return String.format("pointMut(%s: %s:%d %s -> %s)",
                mPointMutation.Gene, mPointMutation.Chromosome,  mPointMutation.Position,
                mPointMutation.Ref, mPointMutation.Alt);
    }

}
