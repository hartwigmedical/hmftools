package com.hartwig.hmftools.svtools.neo;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.fusion.TranscriptUtils.codingBasesToPhase;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionWithin;
import static com.hartwig.hmftools.svtools.neo.NeoEpitopeAnnotator.IM_LOGGER;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.getCodingBases;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.setTranscriptCodingData;
import static com.hartwig.hmftools.svtools.neo.NeoUtils.setTranscriptContext;

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
                    .filter(x -> positionWithin(position(FS_DOWN), x.ExonStart, x.ExonEnd))
                    .findFirst().orElse(null);

            if(exon != null)
            {
                Phases[FS_DOWN] = exon.ExonPhaseEnd;
                ExonRank[FS_DOWN] = exon.ExonRank + 1;
            }
        }
    }

    public void extractCodingBases(final RefGenomeInterface refGenome, int requiredAminoAcids)
    {
        // get the number of bases from the fusion junction as required by the required amino acid count

        // for fusions, the upstream bases include any inserted bases if the region is exonic
        // phasing already takes insert sequence into account, so just adjust for the number of bases required

        int upPhaseOffset = getUpstreamOpenCodonBases();

        int downPhaseOffset = getDownstreamPhaseOffset();

        String codingInsSequence = "";

        if(RegionType[FS_UP] == TranscriptRegionType.EXONIC && !mSvFusion.InsertSequence.isEmpty())
        {
            codingInsSequence = mSvFusion.InsertSequence;
            int upstreamPhase = codingBasesToPhase(Phases[FS_UP] + 1 - codingInsSequence.length());
            upPhaseOffset = getUpstreamOpenCodonBases(upstreamPhase);
        }

        int upRequiredBases = requiredAminoAcids * 3 + upPhaseOffset;

        CodingBases[FS_UP] = getCodingBases(
                refGenome, TransData[FS_UP], chromosome(FS_UP), position(FS_UP), orientation(FS_UP),
                upRequiredBases, true);

        if(!codingInsSequence.isEmpty())
        {
            if(TransData[FS_UP].Strand == POS_STRAND)
                CodingBases[FS_UP] += codingInsSequence;
            else
                CodingBases[FS_UP]= codingInsSequence + CodingBases[FS_UP];
        }

        boolean canStartInExon = RegionType[FS_UP] == TranscriptRegionType.EXONIC;
        int downRequiredBases = requiredAminoAcids * 3 + downPhaseOffset;

        CodingBases[FS_DOWN] = getCodingBases(
                refGenome, TransData[FS_DOWN], chromosome(FS_DOWN), position(FS_DOWN), orientation(FS_DOWN),
                downRequiredBases, canStartInExon);

        IM_LOGGER.trace("ne({}) phased({} up={} down={}) reqBases(up={} down={}) insSeqLen({})",
                toString(), phaseMatched(), upPhaseOffset, downPhaseOffset, upRequiredBases, downRequiredBases, codingInsSequence.length());
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
