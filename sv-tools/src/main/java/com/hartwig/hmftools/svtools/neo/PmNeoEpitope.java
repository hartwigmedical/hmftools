package com.hartwig.hmftools.svtools.neo;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.switchStream;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CodingEffect.MISSENSE;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.getCodingBases;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.setTranscriptCodingData;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.setTranscriptContext;

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
        return (fs == FS_UPSTREAM) == (TransData[fs].Strand == POS_STRAND) ? POS_ORIENT : NEG_ORIENT;
    }

    public int position(int stream)
    {
        // position for a point mutation is defined on the upstream as the base of the mutation
        // for the downstream it is one base before for +ve strand and 1 after for -ve strand
        int pmPosition = mPointMutation.Position;

        if(mIndelBaseDiff >= 0)
        {
            if(stream == FS_UPSTREAM)
                return pmPosition;
            else
                return TransData[FS_UPSTREAM].Strand == POS_STRAND ? pmPosition + 1 : pmPosition - 1;
        }

        // handle the DEL scenario - shift the downstream position by the delete length
        if(stream == FS_UPSTREAM)
        {
            return pmPosition;
        }
        else
        {
            if(TransData[FS_UPSTREAM].Strand == POS_STRAND)
                return pmPosition + abs(mIndelBaseDiff) + 1;
            else
                return pmPosition - 1;
        }
    }

    public String chromosome(int stream)
    {
        return mPointMutation.Chromosome;
    }

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

        TransData[FS_UPSTREAM] = TransData[FS_DOWNSTREAM] = transData;

        int indelBaseDiff = mPointMutation.Alt.length() - mPointMutation.Ref.length();
        int insertedBaseLength = max(indelBaseDiff, 0);

        // set the data for the lower part of the mutation
        int lowerStream = transData.Strand == POS_STRAND ? FS_UPSTREAM : FS_DOWNSTREAM;

        setTranscriptContext(this, transData, position(lowerStream), lowerStream);

        int insSeqLength = lowerStream == FS_UPSTREAM ? insertedBaseLength : 0;

        setTranscriptCodingData(this, transData, position(lowerStream), insSeqLength, lowerStream);

        int upperStream = switchStream(lowerStream);

        // for DELs, set the downstream data as well since it can cross exon-boundaries and/or affect coding bases
        setTranscriptContext(this, transData, position(upperStream), upperStream);

        insSeqLength = upperStream == FS_UPSTREAM ? insertedBaseLength : 0;
        setTranscriptCodingData(this, transData, position(upperStream), insSeqLength, upperStream);
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

        if(TransData[FS_UPSTREAM].Strand == POS_STRAND)
        {
            upPosition = mPointMutation.Position - 1;
            downPosition = mPointMutation.Position + deletedBases + 1;
        }
        else
        {
            upPosition = mPointMutation.Position + deletedBases + 1;
            downPosition = mPointMutation.Position - 1;
        }

        upCodingBases = getCodingBases(
                refGenome, TransData[FS_UPSTREAM], chromosome(FS_UPSTREAM), upPosition, orientation(FS_UPSTREAM),
                upRequiredBases, true);

        // the ref bases was excluded and is now replaced by the alt
        if(TransData[FS_UPSTREAM].Strand == POS_STRAND)
        {
            upCodingBases += mPointMutation.Alt;
        }
        else
        {
            upCodingBases = mPointMutation.Alt + upCodingBases;
        }

        boolean canStartInExon = true; // assumed true for now RegionType[FS_UPSTREAM] == TranscriptRegionType.EXONIC;
        int downRequiredBases = requiredAminoAcids * 3 + downPhaseOffset;

        downCodingBases = getCodingBases(
                refGenome, TransData[FS_DOWNSTREAM], chromosome(FS_DOWNSTREAM), downPosition, orientation(FS_DOWNSTREAM),
                downRequiredBases, canStartInExon);

        CodingBases[FS_UPSTREAM] = upCodingBases;
        CodingBases[FS_DOWNSTREAM] = downCodingBases;
    }

    public String toString()
    {
        return String.format("pointMut(%s: %s:%d %s -> %s)",
                mPointMutation.Gene, mPointMutation.Chromosome,  mPointMutation.Position,
                mPointMutation.Ref, mPointMutation.Alt);
    }

}
