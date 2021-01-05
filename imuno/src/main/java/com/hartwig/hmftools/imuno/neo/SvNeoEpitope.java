package com.hartwig.hmftools.imuno.neo;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_1;
import static com.hartwig.hmftools.common.fusion.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.tickPhaseForward;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.imuno.common.ImunoCommon.IM_LOGGER;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.ALL_TRANS_BASES;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.getDownstreamCodingBases;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.getUpstreamCodingBases;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.setTranscriptCodingData;
import static com.hartwig.hmftools.imuno.neo.NeoUtils.setTranscriptContext;

import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.fusion.TranscriptRegionType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.neo.NeoEpitopeFusion;

public class SvNeoEpitope extends NeoEpitope
{
    private final NeoEpitopeFusion mSvFusion;

    public SvNeoEpitope(final NeoEpitopeFusion fusion)
    {
        super();
        mSvFusion = fusion;
    }

    public int position(int stream)
    {
        return mSvFusion.Positions[stream];
    }

    public String chromosome(int stream)
    {
        return mSvFusion.Chromosomes[stream];
    }
    public String geneName(int stream) { return mSvFusion.GeneNames[stream]; }

    public String variantType()
    {
        return phaseMatched() ? "INFRAME_FUSION" : "OUT_OF_FRAME_FUSION";
    }

    public String variantInfo()
    {
        return String.format("%s:%d:%d-%s:%d:%d",
            mSvFusion.Chromosomes[FS_UP], mSvFusion.Positions[FS_UP], mSvFusion.Orientations[FS_UP],
            mSvFusion.Chromosomes[FS_DOWN], mSvFusion.Positions[FS_DOWN], mSvFusion.Orientations[FS_DOWN]);
    }

    public double copyNumber() { return mSvFusion.JunctionCopyNumber; }

    public boolean phaseMatched()
    {
        if(RegionType[FS_UP] == EXONIC && RegionType[FS_DOWN] == EXONIC)
            return tickPhaseForward(Phases[FS_UP]) == Phases[FS_DOWN];
        else
            return Phases[FS_UP] == Phases[FS_DOWN];
    }

    public void setTranscriptData(final TranscriptData upTransData, final TranscriptData downTransData)
    {
        TransData[FS_UP] = upTransData;
        setTranscriptContext(this, TransData[FS_UP], position(FS_UP), FS_UP);
        int insSeqLength = mSvFusion.InsertSequence.length();
        setTranscriptCodingData(this, TransData[FS_UP], position(FS_UP), insSeqLength, FS_UP);

        TransData[FS_DOWN] = downTransData;

        // if the upstream context is intronic, then skip past any downstream exonic section

        setTranscriptContext(this, TransData[FS_DOWN], position(FS_DOWN), FS_DOWN);
        setTranscriptCodingData(this, TransData[FS_DOWN], position(FS_DOWN), 0, FS_DOWN);

        if(RegionType[FS_UP] == INTRONIC && RegionType[FS_DOWN] == EXONIC)
        {
            final ExonData exon = downTransData.exons().stream()
                    .filter(x -> positionWithin(position(FS_DOWN), x.Start, x.End))
                    .findFirst().orElse(null);

            if(exon != null)
            {
                Phases[FS_DOWN] = exon.PhaseEnd;
                ExonRank[FS_DOWN] = exon.Rank + 1;
            }
        }
    }

    private int getUpstreamOpenCodonBases()
    {
        int phase = Phases[FS_UP];

        if(RegionType[FS_UP] == TranscriptRegionType.EXONIC && !mSvFusion.InsertSequence.isEmpty())
        {
            int insSeqLen = mSvFusion.InsertSequence.length();
            int preInsertPhase = (phase - insSeqLen) % 3;

            if(preInsertPhase >= 0)
                phase = preInsertPhase;
            else if(abs(preInsertPhase) == PHASE_1)
                phase = PHASE_2;
            else
                phase = PHASE_1;
        }

        return getUpstreamOpenCodonBases(phase);
    }

    private int getDownstreamOpenCodonBases()
    {
        int upOpenBases = Phases[FS_UP];
        return upOpenBases == 0 ? 0 : 3 - upOpenBases;
    }

    public void extractCodingBases(final RefGenomeInterface refGenome, int requiredAminoAcids)
    {
        // get the number of bases from the fusion junction as required by the required amino acid count

        // for fusions, the upstream bases include any inserted bases if the region is exonic
        // phasing already takes insert sequence into account, so just adjust for the number of bases required

        int upExtraBases = getUpstreamOpenCodonBases();
        int downExtraBases = getDownstreamOpenCodonBases();

        int upRequiredBases = requiredAminoAcids * 3 + upExtraBases;

        CodingBases[FS_UP] = getUpstreamCodingBases(
                refGenome, TransData[FS_UP], chromosome(FS_UP), position(FS_UP), orientation(FS_UP), upRequiredBases);

        int codingInsSeqLen = 0;

        if(RegionType[FS_UP] == TranscriptRegionType.EXONIC && !mSvFusion.InsertSequence.isEmpty())
        {
            codingInsSeqLen = mSvFusion.InsertSequence.length();

            if(TransData[FS_UP].Strand == POS_STRAND)
                CodingBases[FS_UP] += mSvFusion.InsertSequence;
            else
                CodingBases[FS_UP]= mSvFusion.InsertSequence + CodingBases[FS_UP];
        }

        NovelBaseIndex[FS_UP] = upExtraBases + codingInsSeqLen;
        NovelBaseIndex[FS_DOWN] = downExtraBases;

        boolean canStartInExon = RegionType[FS_UP] == TranscriptRegionType.EXONIC || upExtraBases > 0;
        int downRequiredBases = phaseMatched() ? requiredAminoAcids * 3 + downExtraBases : ALL_TRANS_BASES;

        CodingBases[FS_DOWN] = getDownstreamCodingBases(
                refGenome, TransData[FS_DOWN], chromosome(FS_DOWN), position(FS_DOWN), orientation(FS_DOWN),
                downRequiredBases, canStartInExon, true);

        IM_LOGGER.trace("ne({}) phased({} up={} down={}) reqBases(up={} down={}) insSeqLen({})",
                toString(), phaseMatched(), upExtraBases, downExtraBases, upRequiredBases, downRequiredBases, codingInsSeqLen);
    }

    public String toString()
    {
        return String.format("fusion up(%s: %s:%d:%d) down(%s: %s:%d:%d)",
                mSvFusion.GeneNames[FS_UP], mSvFusion.Chromosomes[FS_UP],
                mSvFusion.Positions[FS_UP], mSvFusion.Orientations[FS_UP],
                mSvFusion.GeneNames[FS_DOWN], mSvFusion.Chromosomes[FS_DOWN],
                mSvFusion.Positions[FS_DOWN], mSvFusion.Orientations[FS_DOWN]);
    }


}
