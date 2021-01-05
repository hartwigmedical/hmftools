package com.hartwig.hmftools.imuno.neo;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.switchStream;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.tickPhaseForward;
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
        // and on the downstream as a first base after the end of the mutated section
        // for negative strand genes, the position is set to the upstream end of the mutated section
        int pmPosition = mPointMutation.Position;
        int modifiedBases = mIndelBaseDiff == 0 ? mPointMutation.Alt.length() : 1;
        int deletedBases = mIndelBaseDiff < 0 ? abs(mIndelBaseDiff) : 0;

        if(posStrand())
        {
            if(stream == FS_UP)
                return pmPosition;
            else
                return pmPosition + deletedBases + modifiedBases;
        }
        else
        {
            // non-INDELs use the start of the mutated bases as the upstream position
            int negPosAdj = mIndelBaseDiff == 0 ? -1 : 0;
            if(stream == FS_UP)
                return pmPosition + deletedBases + modifiedBases + negPosAdj;
            else
                return pmPosition + negPosAdj;
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

    public boolean phaseMatched()
    {
        return (abs(mIndelBaseDiff) % 3) == 0;
    }

    private boolean posStrand() { return TransData[FS_UP].Strand == POS_STRAND; }

    public void setTranscriptData(final TranscriptData upTransData, final TranscriptData downTransData)
    {
        final TranscriptData transData = upTransData;

        TransData[FS_UP] = TransData[FS_DOWN] = transData;

        int insertedBaseLength = max(mIndelBaseDiff, 0);

        // set the data for the lower part of the mutation
        int lowerStream = transData.Strand == POS_STRAND ? FS_UP : FS_DOWN;

        setTranscriptContext(this, transData, position(lowerStream), lowerStream);

        int insSeqLength = 0; // lowerStream == FS_UP ? insertedBaseLength : 0;

        setTranscriptCodingData(this, transData, position(lowerStream), insSeqLength, lowerStream);

        int upperStream = switchStream(lowerStream);

        // for DELs, set the downstream data as well since it can cross exon-boundaries and/or affect coding bases
        setTranscriptContext(this, transData, position(upperStream), upperStream);

        insSeqLength = 0; // upperStream == FS_UP ? insertedBaseLength : 0;
        setTranscriptCodingData(this, transData, position(upperStream), insSeqLength, upperStream);
    }

    private int getUpstreamStartPhase()
    {
        // get the phase at the start of the mutation, taking strand into account
        return Phases[FS_UP];
    }

    private int getUpstreamOpenCodonBases()
    {
        int phase = getUpstreamStartPhase();

        // return the number of open codon bases up to but not including the PM's position
        if(posStrand() || mIndelBaseDiff == 0)
        {
            if(phase == PHASE_1)
                return 0;
            else if(phase == PHASE_2)
                return 1;
            else
                return mIndelBaseDiff == 0 ? 2 : 0; // INDEL leaves the mutation base unchanged, so if it ends a codon, it's not novel
        }
        else
        {
            if(phase == PHASE_1)
                return 1;
            else if(phase == PHASE_2)
                return 2;
            else
                return mIndelBaseDiff == 0 ? 2 : 0;
        }
    }

    private int getDownstreamOpenCodonBases()
    {
        // determine the phase at the base after the mutation - last base of an MNV/SNV, and next base for an INS or DEL
        int mutationTicks;

        if(posStrand())
        {
            mutationTicks = mIndelBaseDiff < 0 ? 1 : mPointMutation.Alt.length();
        }
        else
        {
            // need to go to the phase phase the ALT
            if(mIndelBaseDiff == 0)
                mutationTicks = mPointMutation.Alt.length();
            else if(mIndelBaseDiff > 0)
                mutationTicks = mPointMutation.Alt.length() + 1;
            else
                mutationTicks = 2;
        }

        int postMutationPhase = tickPhaseForward(Phases[FS_UP], mutationTicks);

        if(postMutationPhase == PHASE_0)
            return 1;
        else if(postMutationPhase == PHASE_1)
            return 0;
        else
            return 2;
    }

    public void extractCodingBases(final RefGenomeInterface refGenome, int requiredAminoAcids)
    {
        // get the number of bases from the mutation or junction as required by the required amino acid count
        // the downstream bases start directly after the mutation (ie after any insert, or deleted bases)

        // if the upstream phase is not 0 (ie the end of a codon) then additionally extract enough bases to cover the required amino
        // acids plus the bases of the partial codon - the downstream bases will provide the remainder to complete the codon

        // phasing already takes any inserted bases (ie from an INDEL) into account, so just adjust for the number of bases required
        int upExtraBases = getUpstreamOpenCodonBases();
        int downExtraBases = getDownstreamOpenCodonBases();

        int upRequiredBases = requiredAminoAcids * 3 + upExtraBases;

        String upCodingBases = "";
        String downCodingBases = "";

        int upPosition = position(FS_UP);
        int downPosition = position(FS_DOWN);

        // adjust the position for getting reference bases to make way for the ALT to be added
        if(posStrand())
        {
            upPosition = mPointMutation.Position - 1;
        }
        else
        {
            if(mIndelBaseDiff == 0)
                upPosition += 1; // since would otherwise cross with the SNV/MNV

            downPosition = mPointMutation.Position - 1; // since down pos on -ve strand is the mutation base
        }

        upCodingBases = getUpstreamCodingBases(
                refGenome, TransData[FS_UP], chromosome(FS_UP), upPosition, orientation(FS_UP), upRequiredBases);

        boolean canStartInExon = true; // assumed true for now RegionType[FS_UPSTREAM] == TranscriptRegionType.EXONIC;
        int downRequiredBases = phaseMatched() ? requiredAminoAcids * 3 + downExtraBases : ALL_TRANS_BASES;

        downCodingBases = getDownstreamCodingBases(
                refGenome, TransData[FS_DOWN], chromosome(FS_DOWN), downPosition, orientation(FS_DOWN),
                downRequiredBases, canStartInExon, false);

        // the ref bases was excluded and is now replaced by the alt

        int altLength = mPointMutation.Alt.length();
        int upPhase = getUpstreamStartPhase();

        if(posStrand())
        {
            upCodingBases += mPointMutation.Alt;

            if(mIndelBaseDiff == 0)
            {
                NovelBaseIndex[FS_UP] = upExtraBases + altLength;
            }
            else
            {
                // if the mutation occurs on the last base of a codon, then then novel codon starts immediately after
                if(upPhase == PHASE_0)
                {
                    if(mIndelBaseDiff < 0)
                        NovelBaseIndex[FS_UP] = 0;
                    else
                        NovelBaseIndex[FS_UP] = max(altLength - 1, 0);
                }
                else
                {
                    NovelBaseIndex[FS_UP] = upExtraBases + altLength;
                }
            }

            if(phaseMatched())
                NovelBaseIndex[FS_DOWN] = downExtraBases;
        }
        else
        {
            downCodingBases += mPointMutation.Alt;

            NovelBaseIndex[FS_UP] = upExtraBases;

            if(phaseMatched())
            {
                if(mIndelBaseDiff < 0 && upPhase == PHASE_0)
                    NovelBaseIndex[FS_DOWN] = 0;
                else
                    NovelBaseIndex[FS_DOWN] = downExtraBases + altLength;
            }
        }

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
